// Karl Yerkes
// 2022-01-20

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"
#include "al/graphics/al_Image.hpp"
#include "al/io/al_File.hpp"

#include <cmath>

using namespace al;

#include <fstream>
#include <vector>
using namespace std;

string slurp(string fileName);  // forward declaration

struct AlloApp : App {
  Parameter pointSize{"/pointSize", "", 2.5, 0.1, 3.0};
  Parameter timeStep{"/timeStep", "", 0.1, 0.01, 0.6};
  //

  ShaderProgram pointShader;

  int state = 0;
  double dt_add = 0;

  void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(pointSize);  // add parameter to GUI
    gui.add(timeStep);   // add parameter to GUI
    //
  }

  // a mesh for every style
  Mesh current;
  Mesh original;
  Mesh rgb;
  Mesh hsv;
  Mesh your_style;

  void onCreate() override {
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    current.primitive(Mesh::POINTS);

    auto file = File::currentPath() + "../image-erik-mondrian.png";
    // auto file = File::currentPath() + "../colorful.png";
    auto image = Image(file);
    if (image.width() == 0) {
      cout << "did not load image" << endl;
      exit(1);
    }
    auto aspect_ratio = 1.0f * image.width() / image.height();
    for (int j = 0; j < image.height(); j++) {
      for (int i = 0; i < image.width(); i++) {
        auto pixel = image.at(i, j); // 0-255 (unsigned char / uint8)

        RGB rgb_color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);
        // std::cout << rgb_color.r << " " << rgb_color.g << " " << rgb_color.b << "\n";
        HSV hsv_color = rgb_color;

        current.vertex(1.0 * i / image.width() * aspect_ratio, 1.0 * j / image.height(), 0);
        current.color(rgb_color);
        current.texCoord(0.05, 0);  // s, t

        original.vertex(1.0 * i / image.width() * aspect_ratio, 1.0 * j / image.height(), 0);
        original.color(rgb_color);
        original.texCoord(0.05, 0);  // s, t

        rgb.vertex(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);
        // rgb.vertex(pixel.r / 255.0 - 0.5, pixel.g / 255.0 - 0.5, pixel.b / 255.0 - 0.5);
        rgb.color(rgb_color);
        rgb.texCoord(0.05, 0);  // s, t

        hsv.vertex(hsv_color.s * cos(hsv_color.h * 360 * M_PI / 180), hsv_color.v, hsv_color.s * sin(hsv_color.h * 360 * M_PI / 180));
        hsv.color(rgb_color);
        hsv.texCoord(0.05, 0);  // s, t

        your_style.vertex(tan(pixel.r / 255.0) * pixel.g / 255.0, tan(pixel.g / 255.0) * pixel.b / 255.0, tan(pixel.b / 255.0) * pixel.r / 255.0);
        your_style.color(rgb_color);
        your_style.texCoord(0.05, 0);  // s, t

        // XXX add your RGB and HSV mesh construction here.
      }
    }

    nav().pos(0, 0, 5);
  }

  void onAnimate(double dt) override {
    switch (state) {
      case 1:
        for (int i = 0; i < current.vertices().size(); i++) {
          current.vertices()[i].lerp(original.vertices()[i], dt_add * 0.025);
        }
        dt_add += dt;
        break;
      case 2:
        for (int i = 0; i < current.vertices().size(); i++) {
          current.vertices()[i].lerp(rgb.vertices()[i], dt_add * 0.025);
        }
        dt_add += dt;
        break;
      case 3:
        for (int i = 0; i < current.vertices().size(); i++) {
          current.vertices()[i].lerp(hsv.vertices()[i], dt_add * 0.025);
        }
        dt_add += dt;
        break;
      case 4:
        for (int i = 0; i < current.vertices().size(); i++) {
          current.vertices()[i].lerp(your_style.vertices()[i], dt_add * 0.025);
        }
        dt_add += dt;
        break;
    }

    // XXX accumulate dt to animate transitions between meshes
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == '1') {
      state = 1;
      dt_add = 0;
      // XXX trigger a transition from the current state to the RGB state
    }
    if (k.key() == '2') {
      state = 2;
      dt_add = 0;
    }
    if (k.key() == '3') {
      state = 3;
      dt_add = 0;
    }
    if (k.key() == '4') {
      state = 4;
      dt_add = 0;
    }
    // XXX add more key-based triggers here
    return true;
  }

  void onDraw(Graphics &g) override {
    // XXX no need to change anything in this function
    g.clear(0.3);
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.blending(true);
    g.blendTrans();
    g.depthTesting(true);
    g.draw(current);
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
