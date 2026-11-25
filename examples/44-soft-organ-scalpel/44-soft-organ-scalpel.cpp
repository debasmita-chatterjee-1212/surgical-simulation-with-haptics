//==============================================================================
// soft organs: touch probe and scalpel (shallow surface cuts) with force and deformation readout
//
// All physics lives in surgical_sim.h (tensor-mass model, thesis section 4.4.3);
// organ meshes/materials are in organs.h; the visual instrument geometry is built below.
// This file only (1) reads the mouse/keyboard, (2) calls sim.update(), (3) draws the result
// and the on-screen force / deformation readout.
//
// No haptic device needed: the mouse drives the probe. With a real device, replace
// the "PROBE TARGET FROM MOUSE" block by the device position and send the force
// returned by sim.update() to the device - nothing else changes.
//==============================================================================

// the physics headers must come first (before Windows/GL headers define min/max macros)
#include "organs.h"
#include "scalpel.h"

#include <GL/glew.h>
#include "chai3d.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <cstdio>

using namespace chai3d;
using namespace std;

//------------------------------------------------------------------------------
// SETTINGS
//------------------------------------------------------------------------------
bool fullscreen = false;

// where the probe can move (x = depth towards the organ, y/z = across the screen); set per organ
double workspaceHalfY = 0.13;           // [m]
double workspaceHalfZ = 0.10;           // [m]
double probeXMin = -0.03;
double probeXMax = 0.09;
double probeX = 0.08;                   // start outside the organ; arrow keys push in
double probeKeySpeed = 0.08;            // [m/s] speed of the UP/DOWN arrow keys
double probeScrollStep = 0.004;         // [m] per scroll tick

//------------------------------------------------------------------------------
// VARIABLES
//------------------------------------------------------------------------------
cWorld* world;
cCamera* camera;
cSpotLight* light;
cMesh* tissue;                          // the deformable organ (surface of the tetrahedral mesh)
cMesh* probeMesh;                       // the instrument: a slender probe with a small spherical tip
cMesh *grooveMesh, *slitMesh;           // the cuts: a dark red ribbon (gape) and a black centre line (slit)
cLabel *labForce, *labBar, *labDeform, *labInfo, *labCut, *labKeys;

enum Tool { TOOL_PROBE = 0, TOOL_SCALPEL = 1 };
Tool tool = TOOL_PROBE;
tmass::Scalpel scalpel;

tmass::SurgicalSim sim;
tmass::SimConfig simConfig;
int organId = tmass::ORGAN_LIVER;
double youngBase = 0.0;                 // preset Young's modulus (the [ ] keys scale around it)

double mouseX = 0.0, mouseY = 0.0;
bool showEdges = false;

const double instrumentShaftDiameter = 0.0025; // [m] 2.5 mm

int windowW = 0, windowH = 0;           // framebuffer size (for rendering)
int winW = 1, winH = 1;                 // window size (for mouse mapping)
GLFWwindow* window = nullptr;
bool simulationRunning = false;

//------------------------------------------------------------------------------
// FUNCTIONS
//------------------------------------------------------------------------------
void onErrorCallback(int a_error, const char* a_description);
void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);
void onCursorPosCallback(GLFWwindow* a_window, double a_posX, double a_posY);
void onScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY);
void rebuildTissueMesh();
void updateTissueVertices();
void selectOrgan(int id);
void rebuildProbeMesh();
void setProbeMaterial(double forceFraction);
void rebuildGroove();
void addRibbon(cMesh* m, const std::vector<cVector3d>& P, const std::vector<cVector3d>& N, const std::vector<double>& halfW, double lift);
void updateReadout(double fN, double fMax);

//==============================================================================
int main(int argc, char* argv[])
{
    cout << endl;
    cout << "-----------------------------------" << endl;
    cout << "CHAI3D" << endl;
    cout << "Demo: 40-surgical-sim (soft organs: probe and scalpel - liver, ear, kidney, pancreas)" << endl;
    cout << "-----------------------------------" << endl << endl;
    cout << "Organ:" << endl;
    cout << "  [1] liver   [2] ear cartilage   [3] left kidney   [4] pancreas" << endl;
    cout << "Probe (force and deformation are shown on screen):" << endl;
    cout << "  mouse / touchpad move  - move the probe across the organ" << endl;
    cout << "  UP / DOWN arrow keys   - push the probe in / pull it out (scrolling works too)" << endl;
    cout << "Scalpel:" << endl;
    cout << "  [t] switch probe / scalpel   hold SPACE while pressing in to cut (shallow cuts, at most 4 mm deep)" << endl;
    cout << "  [c] clear all cuts" << endl;
    cout << "Tissue:" << endl;
    cout << "  [ and ]  - softer / stiffer tissue (resets the organ)" << endl;
    cout << "Other keys:" << endl;
    cout << "  [r] reset   [e] mesh edges   [f] full screen   [q] quit" << endl << endl;

    //--------------------------------------------------------------------
    // WINDOW
    //--------------------------------------------------------------------
    glfwSetErrorCallback(onErrorCallback);
    if (!glfwInit())
    {
        cout << "failed to initialize GLFW" << endl;
        return 1;
    }
    glfwWindowHint(GLFW_SAMPLES, 4);
    window = glfwCreateWindow(1000, 700, "CHAI3D - Surgical Sim", nullptr, nullptr);
    if (!window)
    {
        cout << "failed to create window" << endl;
        glfwTerminate();
        return 1;
    }
    glfwSetKeyCallback(window, onKeyCallback);
    glfwSetCursorPosCallback(window, onCursorPosCallback);
    glfwSetScrollCallback(window, onScrollCallback);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // load the real OpenGL function pointers (without this CHAI3D crashes on the first frame)
#ifdef GLEW_VERSION
    if (glewInit() != GLEW_OK)
    {
        cout << "failed to initialize GLEW library" << endl;
        glfwTerminate();
        return 1;
    }
#endif

    //--------------------------------------------------------------------
    // WORLD - CAMERA - LIGHT
    //--------------------------------------------------------------------
    world = new cWorld();
    world->m_backgroundColor.setBlack();

    camera = new cCamera(world);
    world->addChild(camera);
    camera->setClippingPlanes(0.005, 10.0);

    light = new cSpotLight(world);
    camera->addChild(light);
    light->setEnabled(true);
    light->setLocalPos(0.3, 0.3, 0.3);
    light->setDir(-1.0, -1.0, -1.0);

    //--------------------------------------------------------------------
    // TISSUE + PROBE MESHES (created once; selectOrgan() fills them)
    //--------------------------------------------------------------------
    tissue = new cMesh();
    world->addChild(tissue);
    tissue->m_material->m_ambient.set(0.35f, 0.12f, 0.12f, 1.0f);
    tissue->m_material->m_diffuse.set(0.85f, 0.42f, 0.42f, 1.0f);
    tissue->m_material->m_specular.set(0.35f, 0.35f, 0.35f, 1.0f);

    probeMesh = new cMesh();
    world->addChild(probeMesh);
    probeMesh->m_material->m_ambient.set(0.25f, 0.25f, 0.25f, 1.0f);
    probeMesh->m_material->m_diffuse.set(0.75f, 0.75f, 0.80f, 1.0f);
    probeMesh->m_material->m_specular.set(0.9f, 0.9f, 0.9f, 1.0f);

    grooveMesh = new cMesh(); world->addChild(grooveMesh);
    grooveMesh->m_material->m_ambient.set(0.30f, 0.02f, 0.02f, 1.0f);
    grooveMesh->m_material->m_diffuse.set(0.38f, 0.03f, 0.03f, 1.0f);
    grooveMesh->m_material->m_specular.set(0.0f, 0.0f, 0.0f, 1.0f);
    grooveMesh->setUseCulling(false);
    slitMesh = new cMesh(); world->addChild(slitMesh);
    slitMesh->m_material->m_ambient.set(0.02f, 0.0f, 0.0f, 1.0f);
    slitMesh->m_material->m_diffuse.set(0.04f, 0.0f, 0.0f, 1.0f);
    slitMesh->m_material->m_specular.set(0.0f, 0.0f, 0.0f, 1.0f);
    slitMesh->setUseCulling(false);

    // on-screen readout
    labForce  = new cLabel(NEW_CFONT_CALIBRI_28());
    labBar    = new cLabel(NEW_CFONT_CALIBRI_20());
    labDeform = new cLabel(NEW_CFONT_CALIBRI_28());
    labInfo   = new cLabel(NEW_CFONT_CALIBRI_20());
    labCut    = new cLabel(NEW_CFONT_CALIBRI_20());
    labKeys   = new cLabel(NEW_CFONT_CALIBRI_16());
    cLabel* labs[] = { labForce, labBar, labDeform, labInfo, labCut, labKeys };
    for (cLabel* l : labs) { camera->m_frontLayer->addChild(l); l->m_fontColor.setWhite(); }
    labKeys->setText("[1] liver  [2] ear  [3] kidney  [4] pancreas  |  [t] probe / scalpel   SPACE: cut   [c] clear cuts  |  arrows / scroll: depth   [ ]: stiffness   [r] reset   [e] edges   [f] full   [q] quit");

    selectOrgan(tmass::ORGAN_LIVER);

    //--------------------------------------------------------------------
    // MAIN LOOP
    //--------------------------------------------------------------------
    simulationRunning = true;
    double lastTime = glfwGetTime();
    double lastPrint = lastTime;

    while (!glfwWindowShouldClose(window) && simulationRunning)
    {
        const double now = glfwGetTime();
        const double dt = now - lastTime;
        lastTime = now;

        glfwGetFramebufferSize(window, &windowW, &windowH);
        glfwGetWindowSize(window, &winW, &winH);

        // ---- keyboard: arrow keys move the probe in depth ----
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)   probeX -= probeKeySpeed * dt;
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) probeX += probeKeySpeed * dt;
        if (probeX < probeXMin) probeX = probeXMin;
        if (probeX > probeXMax) probeX = probeXMax;

        // ---- PROBE TARGET FROM MOUSE (replace by a haptic device position later) ----
        const double nx = (winW > 0) ? (mouseX / (double)winW) : 0.5;     // 0..1 left to right
        const double ny = (winH > 0) ? (mouseY / (double)winH) : 0.5;     // 0..1 top to bottom
        const tmass::Vec3 target(probeX,
                                 (nx - 0.5) * 2.0 * workspaceHalfY,
                                 (0.5 - ny) * 2.0 * workspaceHalfZ);

        // ---- PHYSICS: returns the force the tissue exerts on the probe ----
        const bool cutKey = (tool == TOOL_SCALPEL) && (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS);
        if (tool == TOOL_SCALPEL) scalpel.beforeStep(sim); else sim.setLoadLimit(-1.0);
        const tmass::Vec3 force = sim.update(dt, target, false);
        if (tool == TOOL_SCALPEL) scalpel.afterStep(sim, cutKey, dt);

        // ---- DRAW ----
        if (sim.consumeSurfaceChanged()) rebuildTissueMesh();
        else                             updateTissueVertices();

        rebuildGroove();

        const tmass::Vec3& p = sim.probe();
        probeMesh->setLocalPos(p.x, p.y, p.z);
        const double fN = force.length();
        const double fMax = (simConfig.contactMaxForce > 0.0) ? simConfig.contactMaxForce : 3.0;
        setProbeMaterial((std::min)(1.0, fN / fMax));
        updateReadout(fN, fMax);

        world->updateShadowMaps(false, false);
        camera->renderView(windowW, windowH);
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (now - lastPrint > 0.25)
        {
            lastPrint = now;
            printf("force %6.3f N   max tissue displacement %6.2f mm   %s   %s      \r",
                   force.length(), sim.maxDisplacement() * 1000.0,
                   sim.inContact() ? "in contact" : "no contact",
                   scalpel.cutting() ? "CUTTING" : "");
            fflush(stdout);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

//------------------------------------------------------------------------------
// switch to another organ: new simulation, camera, workspace and instrument geometry
void selectOrgan(int id)
{
    organId = id;
    simConfig = tmass::organPreset(id);
    youngBase = simConfig.youngModulus;
    sim.init(simConfig);
    scalpel.clear();

    const tmass::OrganView view = tmass::organView(id);
    workspaceHalfY = view.workspaceHalfY;
    workspaceHalfZ = view.workspaceHalfZ;
    probeKeySpeed = view.keySpeed;
    probeScrollStep = view.scrollStep;

    // probe depth range: from well in front of the organ to its rear face
    const double front = sim.boundsMax().x;
    const double extent = (std::max)(sim.boundsMax().y - sim.boundsMin().y, sim.boundsMax().z - sim.boundsMin().z);
    probeXMax = front + 0.5 * extent;
    probeXMin = sim.boundsMin().x;
    probeX = probeXMax;

    camera->set(cVector3d(view.cameraDist, 0.0, 0.12 * view.cameraDist),    // camera position
                cVector3d(0.0, 0.0, 0.0),                                    // look-at
                cVector3d(0.0, 0.0, 1.0));                                   // up

    rebuildProbeMesh();
    rebuildTissueMesh();

    cout << endl << "organ: " << tmass::organName(id) << " - " << sim.activeTets() << " tetrahedra, time step "
         << sim.timeStep() * 1e6 << " us" << endl;
}

//------------------------------------------------------------------------------
void rebuildProbeMesh()
{
    // The instrument enters the organ along -X. Its working tip is at the
    // simulation probe position; the shaft and handle extend away along +X.
    // Geometry is visual only. Contact remains the small spherical proxy in
    // SurgicalSim::applyContact().
    probeMesh->clear();

    if (tool == TOOL_SCALPEL)
    {
        // Scalpel: a flat blade (tip at the origin, pointing along -X, cutting edge on the lower side) and a round handle.
        const double th = 0.0002;                                   // half thickness 0.2 mm
        const double bx[] = { 0.0, 0.0022, 0.0060, 0.0140, 0.0150 };   // blade outline: x ...
        const double by[] = { 0.0, 0.0017, 0.0030, 0.0030, 0.0000 };   // ... and half height (y)
        const int nb = 5;
        std::vector<unsigned int> top, bot;
        for (int i = 0; i < nb; ++i)
        {
            top.push_back(probeMesh->newVertex(cVector3d(bx[i], by[i], +th)));
            bot.push_back(probeMesh->newVertex(cVector3d(bx[i], by[i], -th)));
        }
        for (int i = 1; i + 1 < nb; ++i)                           // two faces (fan around vertex 0 of the outline)
        {
            probeMesh->newTriangle(top[0], top[i], top[i + 1]);
            probeMesh->newTriangle(bot[0], bot[i + 1], bot[i]);
        }
        for (int i = 0; i < nb; ++i)                               // rim
        {
            const int k = (i + 1) % nb;
            probeMesh->newTriangle(top[i], bot[i], bot[k]);
            probeMesh->newTriangle(top[i], bot[k], top[k]);
        }
        // handle: a tapered cylinder along +X
        const double hx0 = 0.0148, hx1 = 0.100, r0 = 0.0016, r1 = 0.0040; const int hs = 14;
        std::vector<unsigned int> ra, rb;
        for (int j = 0; j < hs; ++j)
        {
            const double a = 2.0 * C_PI * j / hs;
            ra.push_back(probeMesh->newVertex(cVector3d(hx0, r0 * cos(a), r0 * sin(a))));
            rb.push_back(probeMesh->newVertex(cVector3d(hx1, r1 * cos(a), r1 * sin(a))));
        }
        for (int j = 0; j < hs; ++j)
        {
            const int k = (j + 1) % hs;
            probeMesh->newTriangle(ra[j], ra[k], rb[k]);
            probeMesh->newTriangle(ra[j], rb[k], rb[j]);
        }
        const unsigned int cap = probeMesh->newVertex(cVector3d(hx1, 0.0, 0.0));
        for (int j = 0; j < hs; ++j) probeMesh->newTriangle(rb[j], rb[(j + 1) % hs], cap);
        probeMesh->computeAllNormals();
        probeMesh->markForUpdate(false);
        setProbeMaterial(0.0);
        return;
    }

    struct ProfilePoint { double x; double radius; };
    const double shaftR = 0.5 * instrumentShaftDiameter;
    std::vector<ProfilePoint> profile;

    {
        // Rounded dissector/probe tip, 2.5 mm shaft, and a modest proximal grip.
        profile.push_back(ProfilePoint{-shaftR, 0.0});
        profile.push_back(ProfilePoint{-0.75 * shaftR, 0.66 * shaftR});
        profile.push_back(ProfilePoint{-0.25 * shaftR, 0.968 * shaftR});
        profile.push_back(ProfilePoint{0.0, shaftR});
        profile.push_back(ProfilePoint{0.001, shaftR});
        profile.push_back(ProfilePoint{0.045, shaftR});
        profile.push_back(ProfilePoint{0.050, 1.5 * shaftR});
        profile.push_back(ProfilePoint{0.052, 2.8 * shaftR});
        profile.push_back(ProfilePoint{0.070, 2.8 * shaftR});
        profile.push_back(ProfilePoint{0.072, 1.5 * shaftR});
        profile.push_back(ProfilePoint{0.073, 0.0});
    }

    const int sides = 16;
    std::vector<std::vector<unsigned int> > rings;
    rings.reserve(profile.size());
    for (size_t i = 0; i < profile.size(); ++i)
    {
        std::vector<unsigned int> ring;
        if (profile[i].radius <= 1e-12)
            ring.push_back(probeMesh->newVertex(cVector3d(profile[i].x, 0.0, 0.0)));
        else
        {
            for (int j = 0; j < sides; ++j)
            {
                const double a = 2.0 * C_PI * j / sides;
                ring.push_back(probeMesh->newVertex(cVector3d(
                    profile[i].x, profile[i].radius * cos(a), profile[i].radius * sin(a))));
            }
        }
        rings.push_back(ring);
    }

    // Join adjacent profile rings. Single-vertex rings create rounded/capped ends.
    for (size_t i = 0; i + 1 < rings.size(); ++i)
    {
        const std::vector<unsigned int>& a = rings[i];
        const std::vector<unsigned int>& b = rings[i + 1];
        if (a.size() == 1 && b.size() == (size_t)sides)
        {
            for (int j = 0; j < sides; ++j)
                probeMesh->newTriangle(a[0], b[(j + 1) % sides], b[j]);
        }
        else if (a.size() == (size_t)sides && b.size() == 1)
        {
            for (int j = 0; j < sides; ++j)
                probeMesh->newTriangle(a[j], a[(j + 1) % sides], b[0]);
        }
        else if (a.size() == (size_t)sides && b.size() == (size_t)sides)
        {
            for (int j = 0; j < sides; ++j)
            {
                const int k = (j + 1) % sides;
                probeMesh->newTriangle(a[j], a[k], b[k]);
                probeMesh->newTriangle(a[j], b[k], b[j]);
            }
        }
    }
    probeMesh->computeAllNormals();
    probeMesh->markForUpdate(false);
    setProbeMaterial(0.0);
}

// the tip shades from steel grey towards red as the force approaches its limit
void setProbeMaterial(double f)
{
    probeMesh->m_material->m_ambient.set(0.25f, 0.25f, 0.25f, 1.0f);
    probeMesh->m_material->m_diffuse.set(float(0.75 + 0.20 * f), float(0.78 - 0.55 * f), float(0.82 - 0.62 * f), 1.0f);
    probeMesh->m_material->m_specular.set(0.9f, 0.9f, 0.9f, 1.0f);
}

// force, force bar, deformation and organ labels (pixel origin bottom-left)
void updateReadout(double fN, double fMax)
{
    char b[160];
    snprintf(b, sizeof b, "Force: %.2f N", fN);
    labForce->setText(b);
    const int n = 30;
    int fill = (int)(n * (std::min)(1.0, fN / fMax) + 0.5);
    std::string bar = "[";
    for (int i = 0; i < n; ++i) bar += (i < fill) ? '|' : '.';
    snprintf(b, sizeof b, "] 0 - %.1f N", fMax);
    labBar->setText(bar + b);
    snprintf(b, sizeof b, "Deformation: %.2f mm (max tissue displacement)", sim.maxDisplacement() * 1000.0);
    labDeform->setText(b);
    labInfo->setText(std::string(tmass::organName(organId)) + (sim.inContact() ? "   -   in contact" : "   -   no contact"));
    if (tool == TOOL_SCALPEL)
    {
        snprintf(b, sizeof b, "Scalpel: %s   cut depth %.1f mm   total cut length %.0f mm%s",
                 scalpel.cutting() ? "CUTTING" : "ready (hold SPACE and press in)", scalpel.depth() * 1000.0,
                 scalpel.length() * 1000.0, scalpel.onCut() ? "   (on an existing cut: weaker tissue)" : "");
        labCut->setText(b);
    }
    else labCut->setText("Probe");
    const int h = windowH;
    labForce->setLocalPos(20, h - 40);
    labBar->setLocalPos(20, h - 66);
    labDeform->setLocalPos(20, h - 104);
    labInfo->setLocalPos(20, h - 134);
    labCut->setLocalPos(20, h - 160);
    labKeys->setLocalPos(20, 14);
}

//------------------------------------------------------------------------------
// the cuts: points anchored to surface triangles, drawn as a ribbon lying just above the (deformed) surface
void addRibbon(cMesh* m, const std::vector<cVector3d>& P, const std::vector<cVector3d>& N, const std::vector<double>& halfW, double lift)
{
    if (P.size() < 2) return;
    unsigned int prevL = 0, prevR = 0;
    for (size_t i = 0; i < P.size(); ++i)
    {
        const cVector3d a = P[(i == 0) ? 0 : i - 1], b = P[(i + 1 == P.size()) ? i : i + 1];
        cVector3d t = b - a;
        cVector3d side = cCross(N[i], t);
        if (side.length() < 1e-12) { if (i > 0) continue; side = cVector3d(0, 1, 0); }
        side.normalize();
        const cVector3d c = P[i] + N[i] * lift;
        const unsigned int l = m->newVertex(c + side * halfW[i]);
        const unsigned int r = m->newVertex(c - side * halfW[i]);
        if (i > 0)
        {
            m->newTriangle(prevL, prevR, r);
            m->newTriangle(prevL, r, l);
        }
        prevL = l; prevR = r;
    }
}

void rebuildGroove()
{
    grooveMesh->clear();
    slitMesh->clear();
    const std::vector<tmass::CutStroke>& strokes = scalpel.strokes();
    if (strokes.empty()) return;
    const std::vector<tmass::Vec3>& D = sim.displayPositions();
    const std::vector<tmass::SurfaceTri>& S = sim.surface();
    for (size_t si = 0; si < strokes.size(); ++si)
    {
        const std::vector<tmass::CutSample>& v = strokes[si].s;
        std::vector<cVector3d> P, N; std::vector<double> w, w2;
        for (size_t k = 0; k < v.size(); ++k)
        {
            const tmass::SurfaceTri& t = S[(size_t)v[k].tri];
            const tmass::Vec3 a = D[(size_t)t.a], b = D[(size_t)t.b], c = D[(size_t)t.c];
            const tmass::Vec3 p = a * v[k].b[0] + b * v[k].b[1] + c * v[k].b[2];
            tmass::Vec3 n = cross(b - a, c - a);
            const double nl = n.length();
            if (!(nl > 1e-14)) continue;
            n = n * (1.0 / nl);                                   // outward
            P.push_back(cVector3d(p.x, p.y, p.z));
            N.push_back(cVector3d(n.x, n.y, n.z));
            w.push_back(0.0006 + 0.35 * v[k].depth);              // the gape widens with the depth of the cut
            w2.push_back(0.00012 + 0.06 * v[k].depth);
        }
        addRibbon(grooveMesh, P, N, w, 0.00030);
        addRibbon(slitMesh, P, N, w2, 0.00045);
    }
    grooveMesh->computeAllNormals(); grooveMesh->markForUpdate(false);
    slitMesh->computeAllNormals(); slitMesh->markForUpdate(false);
}

//------------------------------------------------------------------------------
// rebuild the whole surface mesh (after a cut or a reset)
void rebuildTissueMesh()
{
    tissue->clear();
    const std::vector<tmass::Vec3>& P = sim.displayPositions();
    for (size_t i = 0; i < P.size(); i++) tissue->newVertex(cVector3d(P[i].x, P[i].y, P[i].z));
    const std::vector<tmass::SurfaceTri>& S = sim.surface();
    for (size_t i = 0; i < S.size(); i++)
        tissue->newTriangle((unsigned int)S[i].a, (unsigned int)S[i].b, (unsigned int)S[i].c);
    tissue->computeAllNormals();
    tissue->setShowEdges(showEdges);
    tissue->markForUpdate(false);
}

//------------------------------------------------------------------------------
// move the existing vertices to the current simulated positions
void updateTissueVertices()
{
    const std::vector<tmass::Vec3>& P = sim.displayPositions();
    for (size_t i = 0; i < P.size(); i++)
        tissue->m_vertices->setLocalPos((unsigned int)i, cVector3d(P[i].x, P[i].y, P[i].z));
    tissue->computeAllNormals();
    tissue->markForUpdate(false);
}

//------------------------------------------------------------------------------
void onErrorCallback(int a_error, const char* a_description)
{
    cout << "GLFW Error: " << a_description << endl;
}

//------------------------------------------------------------------------------
void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    if (a_action != GLFW_PRESS) return;

    if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
    {
        simulationRunning = false;
    }
    else if (a_key == GLFW_KEY_1) selectOrgan(tmass::ORGAN_LIVER);
    else if (a_key == GLFW_KEY_2) selectOrgan(tmass::ORGAN_EAR);
    else if (a_key == GLFW_KEY_3) selectOrgan(tmass::ORGAN_KIDNEY);
    else if (a_key == GLFW_KEY_4) selectOrgan(tmass::ORGAN_PANCREAS);
    else if (a_key == GLFW_KEY_R)
    {
        sim.reset();
        scalpel.clear();
        rebuildTissueMesh();
    }
    else if (a_key == GLFW_KEY_T)
    {
        tool = (tool == TOOL_PROBE) ? TOOL_SCALPEL : TOOL_PROBE;
        rebuildProbeMesh();
    }
    else if (a_key == GLFW_KEY_C)
    {
        scalpel.clear();
    }
    else if ((a_key == GLFW_KEY_LEFT_BRACKET) || (a_key == GLFW_KEY_RIGHT_BRACKET))
    {
        // softer ([) or stiffer (]) tissue: re-initialise with the new Young's modulus
        // (limited to 0.1x .. 10x the literature value of the current organ)
        const double penaltyRatio = (simConfig.youngModulus > 0.0) ? simConfig.contactStiffnessPerArea / simConfig.youngModulus : 0.0;
        simConfig.youngModulus *= (a_key == GLFW_KEY_LEFT_BRACKET) ? 0.7 : (1.0 / 0.7);
        if (simConfig.youngModulus < 0.1 * youngBase) simConfig.youngModulus = 0.1 * youngBase;
        if (simConfig.youngModulus > 10.0 * youngBase) simConfig.youngModulus = 10.0 * youngBase;
        simConfig.contactStiffnessPerArea = penaltyRatio * simConfig.youngModulus;
        sim.init(simConfig);
        scalpel.clear();
        rebuildTissueMesh();
        cout << endl << "Young's modulus = " << simConfig.youngModulus * 1e-3 << " kPa, time step "
             << sim.timeStep() * 1e6 << " us" << endl;
    }
    else if (a_key == GLFW_KEY_E)
    {
        showEdges = !showEdges;
        tissue->setShowEdges(showEdges);
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
        else
        {
            glfwSetWindowMonitor(a_window, nullptr, 50, 50, 1000, 700, 0);
        }
    }
}

//------------------------------------------------------------------------------
void onCursorPosCallback(GLFWwindow* a_window, double a_posX, double a_posY)
{
    mouseX = a_posX;
    mouseY = a_posY;
}

//------------------------------------------------------------------------------
void onScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY)
{
    probeX -= a_offsetY * probeScrollStep;             // scroll forward = push in (towards -x)
    if (probeX < probeXMin) probeX = probeXMin;
    if (probeX > probeXMax) probeX = probeXMax;
}
