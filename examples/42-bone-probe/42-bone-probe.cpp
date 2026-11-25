//==============================================================================
// 42-bone-probe.cpp   (rigid temporal bone, touch probe: contact force and deformation only)
//
// A small sphere (the instrument) is moved over a rigid bone model with the mouse. The bone is NOT cut.
// The contact force and the contact deformation are shown continuously on screen and written to
// bone_probe_log.csv.
//
// Physics is in voxel_bone.h (class vbone::VoxelBone, probe-only mode):
//   - the sphere cannot enter the bone; the hand (mouse target) can, through a virtual spring (800 N/m),
//   - the force felt is that spring force along the surface normal, limited to 3.3 N (haptic device peak),
//   - the deformation is the Hertz indentation of the sphere into the bone (+ whole-bone compliance, when available).
//
// Bone model: BodyParts3D "Left temporal bone" (FJ580, FMA52739),
//   (c) The Database Center for Life Science, licensed under CC Attribution 4.0 International.
// It is a generic solid surface: no air cells, no nerve, no sinus, no ossicles.
//
// No haptic device needed. With a real device, replace the "PROBE TARGET FROM MOUSE" block by the device
// position and send the returned force to the device.
//==============================================================================

// physics header first (before Windows/GL headers define min/max macros)
#include "voxel_bone.h"

#include <GL/glew.h>
#include "chai3d.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

using namespace chai3d;
using namespace std;

//------------------------------------------------------------------------------
// SETTINGS
//------------------------------------------------------------------------------
const char* kVoxelFile   = "temporal_bone_voxels.bin";
const char* kSurfaceFile = "temporal_bone_surface.bin";

double tipRadius   = 0.0015;            // [m] instrument sphere (keys [ and ] change it)
const double tipRadiusMin = 0.0005, tipRadiusMax = 0.0030;
double probeZ      = 0.075;             // [m] height of the hand above the grid floor; arrow keys / scroll change it
const double probeZMin = 0.000, probeZMax = 0.080;
const double probeKeySpeed = 0.05;      // [m/s]

//------------------------------------------------------------------------------
// VARIABLES
//------------------------------------------------------------------------------
cWorld* world;
cCamera* camera;
cSpotLight* light;
cMesh* boneMesh;
cMesh* tipMesh;
cLabel *labForce, *labDeform, *labInfo, *labBar, *labKeys, *labCredit;

vbone::VoxelBone bone;
vbone::BurrConfig cfg;
double boxX = 0.1, boxY = 0.08, boxZ = 0.06;        // grid size [m], filled when the data is loaded

double mouseX = 0.0, mouseY = 0.0;
int windowW = 0, windowH = 0, winW = 1, winH = 1;
GLFWwindow* window = nullptr;
bool simulationRunning = false, fullscreen = false, showEdges = false, tiltedView = false;

//------------------------------------------------------------------------------
void onErrorCallback(int a_error, const char* a_description) { cout << "GLFW Error: " << a_description << endl; }
void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);
void onCursorPosCallback(GLFWwindow* a_window, double a_posX, double a_posY) { mouseX = a_posX; mouseY = a_posY; }
void onScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY);
void setCamera();
void rebuildTip();

// look for a data file next to the executable, in the working folder, next to the source file, or in argv[1]
std::string findFile(const char* name, const std::string& hint)
{
    std::vector<std::string> dirs;
    if (!hint.empty()) dirs.push_back(hint + "/");
    dirs.push_back("");
    dirs.push_back("../../examples/GLFW/42-bone-probe/");
    dirs.push_back("../../../examples/GLFW/42-bone-probe/");
    std::string src = __FILE__;
    size_t k = src.find_last_of("/\\");
    if (k != std::string::npos) dirs.push_back(src.substr(0, k + 1));
    for (size_t i = 0; i < dirs.size(); i++)
    {
        std::string p = dirs[i] + name;
        std::ifstream f(p.c_str(), std::ios::binary);
        if (f.good()) return p;
    }
    return "";
}

//==============================================================================
int main(int argc, char* argv[])
{
    cout << endl << "-----------------------------------" << endl;
    cout << "CHAI3D" << endl << "Demo: 42-bone-probe (rigid bone, contact force and deformation)" << endl;
    cout << "-----------------------------------" << endl << endl;
    cout << "  mouse / touchpad move  - move the sphere across the bone" << endl;
    cout << "  UP / DOWN arrows (or scroll) - push the hand down / lift it" << endl;
    cout << "  [ and ]  - smaller / larger sphere" << endl;
    cout << "  [v] top / tilted view   [e] mesh edges   [r] reset   [f] full screen   [q] quit" << endl << endl;

    //--------------------------------------------------------------------
    // DATA
    //--------------------------------------------------------------------
    const std::string hint = (argc > 1) ? argv[1] : "";
    const std::string voxPath = findFile(kVoxelFile, hint), surfPath = findFile(kSurfaceFile, hint);
    if (voxPath.empty() || surfPath.empty())
    {
        cout << "cannot find " << kVoxelFile << " / " << kSurfaceFile << " - copy them next to the exe, or pass the folder as the first argument." << endl;
        return 1;
    }
    int nx = 0, ny = 0, nz = 0; double hmm = 0.0;
    std::vector<float> rho;
    {
        std::ifstream f(voxPath.c_str(), std::ios::binary);
        f.read((char*)&nx, 4); f.read((char*)&ny, 4); f.read((char*)&nz, 4); f.read((char*)&hmm, 8);
        std::vector<uint8_t> raw((size_t)nx * ny * nz);
        f.read((char*)raw.data(), (std::streamsize)raw.size());
        if (!f) { cout << "voxel file is damaged: " << voxPath << endl; return 1; }
        rho.resize(raw.size());
        for (size_t i = 0; i < raw.size(); i++) rho[i] = raw[i] * (1.0f / 255.0f);
    }
    const double h = hmm * 1e-3;
    boxX = nx * h; boxY = ny * h; boxZ = nz * h;

    cfg.probeOnly = true;
    cfg.burrRadius = tipRadius;
    cfg.logPath = "bone_probe_log.csv";
    bone.init(nx, ny, nz, h, vbone::V3(0, 0, 0), rho, cfg);
    cout << "bone: " << nx << " x " << ny << " x " << nz << " voxels of " << hmm << " mm, solid volume "
         << bone.solidVolumeMM3() * 1e-3 << " cm3" << endl;

    //--------------------------------------------------------------------
    // WINDOW
    //--------------------------------------------------------------------
    glfwSetErrorCallback(onErrorCallback);
    if (!glfwInit()) { cout << "failed to initialize GLFW" << endl; return 1; }
    glfwWindowHint(GLFW_SAMPLES, 4);
    window = glfwCreateWindow(1100, 760, "CHAI3D - Bone probe (rigid temporal bone)", nullptr, nullptr);
    if (!window) { cout << "failed to create window" << endl; glfwTerminate(); return 1; }
    glfwSetKeyCallback(window, onKeyCallback);
    glfwSetCursorPosCallback(window, onCursorPosCallback);
    glfwSetScrollCallback(window, onScrollCallback);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
#ifdef GLEW_VERSION
    if (glewInit() != GLEW_OK) { cout << "failed to initialize GLEW library" << endl; glfwTerminate(); return 1; }
#endif

    //--------------------------------------------------------------------
    // WORLD - CAMERA - LIGHT
    //--------------------------------------------------------------------
    world = new cWorld();
    world->m_backgroundColor.set(0.08f, 0.09f, 0.11f);
    camera = new cCamera(world);
    world->addChild(camera);
    setCamera();
    camera->setClippingPlanes(0.01, 10.0);

    light = new cSpotLight(world);
    camera->addChild(light);
    light->setEnabled(true);
    light->setLocalPos(0.3, 0.3, 0.3);
    light->setDir(-1.0, -1.0, -1.0);

    //--------------------------------------------------------------------
    // BONE SURFACE (smooth marching-cubes mesh made by make_bone_voxels.py)
    //--------------------------------------------------------------------
    boneMesh = new cMesh();
    world->addChild(boneMesh);
    {
        std::ifstream f(surfPath.c_str(), std::ios::binary);
        int32_t nv = 0, nf = 0;
        f.read((char*)&nv, 4); f.read((char*)&nf, 4);
        std::vector<float> V((size_t)nv * 3);
        std::vector<uint32_t> T((size_t)nf * 3);
        f.read((char*)V.data(), (std::streamsize)(V.size() * 4));
        f.read((char*)T.data(), (std::streamsize)(T.size() * 4));
        if (!f) { cout << "surface file is damaged: " << surfPath << endl; return 1; }
        for (int i = 0; i < nv; i++) boneMesh->newVertex(cVector3d(V[3 * i], V[3 * i + 1], V[3 * i + 2]));
        for (int i = 0; i < nf; i++) boneMesh->newTriangle(T[3 * i], T[3 * i + 1], T[3 * i + 2]);
        cout << "display surface: " << nv << " vertices, " << nf << " triangles" << endl;
    }
    boneMesh->computeAllNormals();
    boneMesh->m_material->m_ambient.set(0.42f, 0.38f, 0.30f, 1.0f);
    boneMesh->m_material->m_diffuse.set(0.93f, 0.88f, 0.76f, 1.0f);
    boneMesh->m_material->m_specular.set(0.20f, 0.20f, 0.18f, 1.0f);
    boneMesh->m_material->setShininess(20);
    boneMesh->setUseDisplayList(true);

    tipMesh = new cMesh();
    world->addChild(tipMesh);
    rebuildTip();

    //--------------------------------------------------------------------
    // ON-SCREEN READOUT
    //--------------------------------------------------------------------
    cFontPtr font = NEW_CFONT_CALIBRI_20();
    cFontPtr fontSmall = NEW_CFONT_CALIBRI_16();
    labForce  = new cLabel(font);      camera->m_frontLayer->addChild(labForce);
    labBar    = new cLabel(font);      camera->m_frontLayer->addChild(labBar);
    labDeform = new cLabel(font);      camera->m_frontLayer->addChild(labDeform);
    labInfo   = new cLabel(fontSmall); camera->m_frontLayer->addChild(labInfo);
    labKeys   = new cLabel(fontSmall); camera->m_frontLayer->addChild(labKeys);
    labCredit = new cLabel(fontSmall); camera->m_frontLayer->addChild(labCredit);
    labForce->m_fontColor.setWhite(); labBar->m_fontColor.setWhite(); labDeform->m_fontColor.setWhite();
    labInfo->m_fontColor.setGrayLight(); labKeys->m_fontColor.setGrayLight(); labCredit->m_fontColor.setGrayLight();
    labKeys->setText("mouse: move   UP/DOWN: push / lift   [ ]: tip size   v: view   r: reset   q: quit");
    labCredit->setText("Bone: BodyParts3D FJ580 (c) DBCLS, CC BY 4.0. Rigid, uncut. E = 10 GPa and nu = 0.3 are placeholders.");

    //--------------------------------------------------------------------
    // MAIN LOOP
    //--------------------------------------------------------------------
    simulationRunning = true;
    double lastTime = glfwGetTime(), lastPrint = lastTime;
    vbone::V3 force;

    while (!glfwWindowShouldClose(window) && simulationRunning)
    {
        const double now = glfwGetTime();
        double dt = now - lastTime;
        lastTime = now;
        if (dt > 0.05) dt = 0.05;

        glfwGetFramebufferSize(window, &windowW, &windowH);
        glfwGetWindowSize(window, &winW, &winH);

        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)   probeZ -= probeKeySpeed * dt;
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) probeZ += probeKeySpeed * dt;
        if (probeZ < probeZMin) probeZ = probeZMin;
        if (probeZ > probeZMax) probeZ = probeZMax;

        // ---- PROBE TARGET FROM MOUSE (replace by a haptic device position later) ----
        const double nxm = (winW > 0) ? (mouseX / (double)winW) : 0.5;     // 0..1 left to right
        const double nym = (winH > 0) ? (mouseY / (double)winH) : 0.5;     // 0..1 top to bottom
        const double margin = 0.004;
        const vbone::V3 target(margin + nxm * (boxX - 2 * margin),
                               margin + (1.0 - nym) * (boxY - 2 * margin),
                               probeZ);

        // ---- PHYSICS: returns the force the bone exerts on the instrument ----
        bone.config().burrRadius = tipRadius;
        force = bone.update(dt, target);

        // ---- DRAW ----
        const vbone::V3& p = bone.burr();
        tipMesh->setLocalPos(p.x, p.y, p.z);
        const double fr = (std::min)(1.0, force.length() / cfg.maxForce);   // 0..1
        tipMesh->m_material->m_diffuse.set((float)(0.75 + 0.25 * fr), (float)(0.75 - 0.55 * fr), (float)(0.80 - 0.60 * fr), 1.0f);

        const double F = force.length();
        char buf[256];
        snprintf(buf, sizeof(buf), "Contact force: %5.2f N   (device limit %.1f N, spring force %.2f N)", F, cfg.maxForce, bone.rawForceN());
        labForce->setText(buf);
        int nbar = (int)(F / cfg.maxForce * 30.0 + 0.5);
        std::string bar = "[";
        for (int i = 0; i < 30; i++) bar += (i < nbar) ? '#' : '-';
        bar += "]";
        labBar->setText(bar);
        snprintf(buf, sizeof(buf), "Contact deformation: %6.2f um  (Hertz indentation %.2f um)", bone.contactDeformationUM(), bone.hertzIndentationUM());
        labDeform->setText(buf);
        snprintf(buf, sizeof(buf), "tip radius %.1f mm   hand beyond surface %.2f mm   %s", tipRadius * 1e3, bone.lastDepthMM(),
                 bone.inContact() ? "IN CONTACT" : "free");
        labInfo->setText(buf);
        labForce->setLocalPos(20, windowH - 40);
        labBar->setLocalPos(20, windowH - 66);
        labDeform->setLocalPos(20, windowH - 92);
        labInfo->setLocalPos(20, windowH - 116);
        labKeys->setLocalPos(20, 40);
        labCredit->setLocalPos(20, 16);

        world->updateShadowMaps(false, false);
        camera->renderView(windowW, windowH);
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (now - lastPrint > 0.25)
        {
            lastPrint = now;
            printf("force %5.2f N   deformation %6.2f um   %s      \r", F, bone.contactDeformationUM(), bone.inContact() ? "contact" : "free   ");
            fflush(stdout);
        }
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

//------------------------------------------------------------------------------
void setCamera()
{
    const double cx = 0.5 * boxX, cy = 0.5 * boxY, cz = 0.5 * boxZ;
    if (!tiltedView)
        camera->set(cVector3d(cx, cy, cz + 0.22), cVector3d(cx, cy, cz), cVector3d(0.0, 1.0, 0.0));
    else
        camera->set(cVector3d(cx, cy - 0.15, cz + 0.17), cVector3d(cx, cy, cz), cVector3d(0.0, 0.0, 1.0));
}

//------------------------------------------------------------------------------
void rebuildTip()
{
    tipMesh->clear();
    cCreateSphere(tipMesh, tipRadius);
    tipMesh->m_material->m_ambient.set(0.25f, 0.25f, 0.25f, 1.0f);
    tipMesh->m_material->m_diffuse.set(0.75f, 0.75f, 0.80f, 1.0f);
    tipMesh->m_material->m_specular.set(0.9f, 0.9f, 0.9f, 1.0f);
}

//------------------------------------------------------------------------------
void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    if (a_action != GLFW_PRESS) return;
    if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q)) simulationRunning = false;
    else if (a_key == GLFW_KEY_R) { bone.reset(); probeZ = probeZMax; }
    else if (a_key == GLFW_KEY_V) { tiltedView = !tiltedView; setCamera(); }
    else if (a_key == GLFW_KEY_E) { showEdges = !showEdges; boneMesh->setShowEdges(showEdges); }
    else if ((a_key == GLFW_KEY_LEFT_BRACKET) || (a_key == GLFW_KEY_RIGHT_BRACKET))
    {
        tipRadius *= (a_key == GLFW_KEY_LEFT_BRACKET) ? 0.8 : 1.25;
        if (tipRadius < tipRadiusMin) tipRadius = tipRadiusMin;
        if (tipRadius > tipRadiusMax) tipRadius = tipRadiusMax;
        rebuildTip();
    }
    else if (a_key == GLFW_KEY_F)
    {
        fullscreen = !fullscreen;
        if (fullscreen)
        {
            GLFWmonitor* monitor = glfwGetPrimaryMonitor();
            const GLFWvidmode* mode = glfwGetVideoMode(monitor);
            glfwSetWindowMonitor(a_window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
        }
        else glfwSetWindowMonitor(a_window, nullptr, 50, 50, 1100, 760, 0);
    }
}

//------------------------------------------------------------------------------
void onScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY)
{
    probeZ -= a_offsetY * 0.003;
    if (probeZ < probeZMin) probeZ = probeZMin;
    if (probeZ > probeZMax) probeZ = probeZMax;
}
