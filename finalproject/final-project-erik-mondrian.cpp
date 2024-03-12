// Karl Yerkes
// 2022-01-20

#include "al/app/al_DistributedApp.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"
#include <cmath>
#include <math.h>

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}

string slurp(string fileName);  // forward declaration

struct CommonState {
  // char data[60000]; // larger than a UDP packet
  // int frame;
  // float signal;
  float pointSize;
  float timeStep;
  float dragFactor;
  float radius;
  float spring_stiffness;
  float charge;
  bool freeze;
  bool love;
  bool rand_force;
};

struct AlloApp : DistributedAppWithState<CommonState> {
  Parameter pointSize{"/pointSize", "", 1.0, 0.0, 2.0};
  Parameter timeStep{"/timeStep", "", 0.3, 0.01, 0.6};
  Parameter dragFactor{"/dragFactor", "", 0.45, 0.0, 0.9};
  Parameter radius{"/radius", "", 3.0, 0.5, 10.0};
  Parameter spring_stiffness{"/springStiffness", "", 0.75, 0.1, 1.0};
  Parameter charge{"/chargeConstant", "", 0.01, 0.01, 1.0};
  //

  ShaderProgram pointShader;

  //  simulation state
  Mesh mesh;  // position *is inside the mesh* mesh.vertices() are the positions
  vector<Vec3f> velocity;
  vector<Vec3f> force;
  vector<float> mass;

  void onInit() override {
    // set up GUI
    auto cuttleboneDomain =
    CuttleboneStateSimulationDomain<CommonState>::enableCuttlebone(this);
    if (!cuttleboneDomain) {
      std::cerr << "ERROR: Could not start Cuttlebone. Quitting." << std::endl;
      quit();
    }

    if (isPrimary()) {
      auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
      auto &gui = GUIdomain->newGUI();
      gui.add(pointSize);  // add parameter to GUI
      gui.add(timeStep);   // add parameter to GUI
      gui.add(dragFactor);   // add parameter to GUI
      gui.add(radius);
      gui.add(spring_stiffness);
      gui.add(charge);
    }

    if (isPrimary()) {
      state().pointSize = pointSize;
      state().timeStep = timeStep;
      state().dragFactor = dragFactor;
      state().radius = radius;
      state().spring_stiffness = spring_stiffness;
      state().charge = charge;
      state().freeze = freeze;
      state().love = love;
      state().rand_force = rand_force;
    }
    //
  }

  void onCreate() override {
    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // set initial conditions of the simulation
    //

    // c++11 "lambda" function
    auto randomColor = []() { return HSV(rnd::uniform(), 1.0f, 1.0f); };

    mesh.primitive(Mesh::POINTS);
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?
    for (int _ = 0; _ < 1000; _++) {
      mesh.vertex(randomVec3f(5));
      mesh.color(randomColor());

      // float m = rnd::uniform(3.0, 0.5);
      float m = 3 + rnd::normal() / 2;
      if (m < 0.5) m = 0.5;
      mass.push_back(m);

      // using a simplified volume/size relationship
      mesh.texCoord(pow(m, 1.0f / 3), 0);  // s, t

      // separate state arrays
      velocity.push_back(randomVec3f(0.1));
      force.push_back(randomVec3f(1));
    }

    // nav().pos(0, 0, 15);
    nav().pos(0, 0, 0);
  }

  bool freeze = false;
  bool love = false;
  bool rand_force = false;

  void onAnimate(double dt) override {
    if (isPrimary()) {
      state().pointSize = pointSize;
      state().timeStep = timeStep;
      state().dragFactor = dragFactor;
      state().radius = radius;
      state().spring_stiffness = spring_stiffness;
      state().charge = charge;
      state().freeze = freeze;
      state().love = love;
      state().rand_force = rand_force;
    }

    if (state().freeze) return;

    // Calculate forces

    // cout << "X=" << force[0].x << ", Y=" << force[0].y << ", Z=" << force[0].z << "\n";

    // XXX you put code here that calculates gravitational forces and sets
    // accelerations These are pair-wise. Each unique pairing of two particles
    // These are equal but opposite: A exerts a force on B while B exerts that
    // same amount of force on A (but in the opposite direction!) Use a nested
    // for loop to visit each pair once The time complexity is O(n*n)
    //
    // Vec3f has lots of operations you might use...
    // • +=
    // • -=
    // • +
    // • -
    // • .normalize() ~ Vec3f points in the direction as it did, but has length 1
    // • .normalize(float scale) ~ same but length `scale`
    // • .mag() ~ length of the Vec3f
    // • .magSqr() ~ squared length of the Vec3f
    // • .dot(Vec3f f) 
    // • .cross(Vec3f f)

    // drag
    for (int i = 0; i < velocity.size(); i++) {
      force[i] += - velocity[i] * state().dragFactor;
    }

    if (state().rand_force)
      for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        force[i] += randomVec3f(1);
      }

    // Integration
    //
    vector<Vec3f> &position(mesh.vertices());

    float current_radius;
    for (int i = 0; i < velocity.size(); i++) {
      current_radius = sqrt(pow(position[i].x, 2) + pow(position[i].y, 2) + pow(position[i].z, 2));
      force[i].x += -state().spring_stiffness * (current_radius - state().radius) * position[i].x / current_radius;
      force[i].y += -state().spring_stiffness * (current_radius - state().radius) * position[i].y / current_radius;
      force[i].z += -state().spring_stiffness * (current_radius - state().radius) * position[i].z / current_radius;
    }

    // Vec3f distance_vector;
    for (int i = 0; i < velocity.size(); i++)
      for (int j = i + 1; j < velocity.size(); j++) {
        // if (j == i)
        //   continue;
        // distance_vector = (position[i] - position[j]).normalize();
        if (state().love && i % 10 == 0) {
          force[i] += (position[i] - position[j]).normalize(state().charge * -1 / (position[i] - position[j]).magSqr());
        }
        else {
          force[i] += (position[i] - position[j]).normalize(state().charge * 1 / (position[i] - position[j]).magSqr());
        }
        force[j] += (position[j] - position[i]).normalize(state().charge * 1 / (position[j] - position[i]).magSqr());
        // force[i] += (charge * charge / (4 * M_PI * 1)) * ((position[i] - position[j]).mag() / pow((position[i] - position[j]).mag(), 3));
      }

    for (int i = 0; i < velocity.size(); i++) {
      // "semi-implicit" Euler integration
      velocity[i] += force[i] / mass[i] * state().timeStep;
      position[i] += velocity[i] * state().timeStep;
    }

    // clear all accelerations (IMPORTANT!!)
    for (auto &a : force) a.set(0);
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '2') {
      love = !love;
    }

    if (k.key() == '1') {
      rand_force = true;
      // introduce some "random" forces
      /* for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        force[i] += randomVec3f(1);
      } */
    }

    return true;
  }

  bool onKeyUp(const Keyboard &k) override {
    if (k.key() == '1')
      rand_force = false;

    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear(0.3);
    g.shader(pointShader);
    g.shader().uniform("pointSize", state().pointSize / 100);
    g.blending(true);
    g.blendTrans();
    g.depthTesting(true);
    g.draw(mesh);
  }
};

int main() {
  AlloApp app;
  app.configureAudio(48000, 512, 2, 0);
  app.start();
}

string slurp(string fileName) {
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}