#include "points3d.hpp"
#include <matplotlib-cpp/matplotlibcpp.h>
namespace plt = matplotlibcpp;

class afm{
  private:
  points3d surface;
  points3d probe;
  float max_surface = -1000000;
  float min_probe = 1000000;
  float LennardJonesForce(vector<vector<float>> p, vector<vector<float>> s, float dx, float dy, float z);
  public:
  afm(points3d s, points3d p) :surface(s),probe(p){
    // probeの下限の取得
    for(auto &p_point : probe.points){
      if(min_probe > p_point[2]){
         min_probe = p_point[2];
      }
    }
    // surface の上限の取得
    for(auto &s_point : surface.points){
      if(max_surface < s_point[2]){
        max_surface = s_point[2];
      }
    }
  };
  // 画像出力系
  void scanXY(float z);
  void scanXYZ();
};

// 制御なしのconstant hightでの計測
void afm::scanXY(float z){
  //100*100の画像を返せば良さげか?
  int nx = 100;
  int ny = 100;
  float dx = 0.1;
  float dy = 0.1;
  vector<vector<float>> img(nx,vector<float>(ny));
  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      img[x][y] = LennardJonesForce(probe.points, surface.points, x*dx, y*dy, z);
    }
  }
  //出力
  vector<float> plot_img(nx*ny);
  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      plot_img.at(ny * y + x) = img[x][y];
    }
  }
  const float* zptr = &(plot_img[0]);
  const int colors = 1;
  plt::imshow(zptr, nx, ny, colors);
  plt::show();
};
//表面から1~6の範囲でスキャン
void afm::scanXYZ(){
  int nx = 100;
  int ny = 100;
  int nz = 30;
  float dx = 0.1;
  float dy = 0.1;
  float dz = 0.1;
  //surface原子の位置;
  vector<float> px,py,pz;
  
  plt::scatter(px,py,pz);
  plt::save("./img/surface.png");

  for(int z=20;z<nz+20;z++){
    for(auto point : surface.points){
      cout<<point[0]<<","<<point[1]<<","<<point[2]<<endl;
    }
    for(int y=0;y<ny;y++){
      for(int x=0;x<nx;x++){
        float f = LennardJonesForce(probe.points, surface.points, x*dx, y*dy, z*dz);
        cout<<x<<","<<y<<","<<z<<","<<f<<endl;
      }
    }
  }
};

float afm::LennardJonesForce(vector<vector<float>> probe, vector<vector<float>> surface, float dx, float dy, float z){
  float F = 0;
  float eV = 1.609e-19;
  float k = 1.38062e-23;
  float m0 = 1e-26;
  float d = 1e-10;
  float sigma = 2.951;
  float epsilon = 0.2294;
  for(auto p: probe){
    // プローブ先の移動を擬似的に表現
    p[0] += dx;
    p[1] += dy;
    // 高さ制御
    p[2] = p[2]-min_probe+z;
    for(auto s: surface){
      float distance = pow(pow(p[0] - s[0], 2) + pow(p[1] - s[1], 2) + pow(p[2] - s[2], 2), 0.5);// 純粋な距離
      float distance_z = p[2] - s[2];
      // LennardJonesの式に入れて計算
      //F = -48*epsilon*((sigma/r)**12 - 0.5 * (sigma/r)**6)*R[2]/(r**2)
      F += -48 * epsilon * (pow(sigma / distance,12) - 0.5 * pow(sigma / distance, 6)) * distance_z / pow(distance, 2);
    }
  }
  return F;
}


int main(){
  points3d surface(0, true, 4.078);
  points3d probe(1, true, 4.078);
  //平面の生成
  surface.face_centered_cubic(4.1,4.1,4.1);
  surface.trim([=](tuple<float,float,float> xyz){
    float x = get<0>(xyz);
    float y = get<1>(xyz);
    float z = get<2>(xyz);
    return x<=3 && y<=3 && z<=3;
  });
  //プローブの生成
  const float probe_r = 10; //曲率
  probe.face_centered_cubic(2*probe_r+10,2*probe_r+10,2*probe_r+10); 
  probe.shift(-probe_r,-probe_r,-probe_r);
  //平面から2nm離れた原子の影響力はほぼないと考えて、半球状に切る
  probe.trim([=](tuple<float,float,float> xyz){
    float x = get<0>(xyz);
    float y = get<1>(xyz);
    float z = get<2>(xyz);
    float r = probe_r;
    return r*r > x*x + y*y + z*z && z <= 0;
  });
  
  //原子間力の計算を行う
  afm AFM(surface, probe);
  // //制御なしのXY測定
  //AFM.scanXY(1);
  AFM.scanXYZ();
  // //fm変調での取得画像
  // AFM.modulationXY(z):
}