import processing.serial.*;
import ddf.minim.*;

/*
Arena is 1000 x 500 x 500
Center of the arena is (0, 0, 0)


                +Z   +Y
                 |  /
                 | /
                 |/
 ----------------+------------------ +X
                /|
               / |
              /  |
              
Each player has a constant x
and can change their y and z coordinats

*/

class Coordinat
{
  Coordinat(float x, float y, float z)
  {
    m_x = x;
    m_y = y;
    m_z = z;
  }
  float GetX() { return m_x; }
  float GetY() { return m_y; }
  float GetZ() { return m_z; }
 
  void SetX(float x) { m_x = x; }
  void SetY(float y) { m_y = y; }
  void SetZ(float z) { m_z = z; }

  float m_x;
  float m_y;
  float m_z;
};

class Vector
{
  Vector(float dx, float dy, float dz)
  {
    m_dx = dx;
    m_dy = dy;
    m_dz = dz;
  }
  float GetDx() { return m_dx; }
  float GetDy() { return m_dy; }
  float GetDz() { return m_dz; }
  void SetDx(float dx) { m_dx = dx; }
  void SetDy(float dy) { m_dy = dy; }
  void SetDz(float dz) { m_dz = dz; }
  float m_dx;
  float m_dy;
  float m_dz;
};

class Arena
{
  Arena(float minx, float maxx, float miny, float maxy, float minz, float maxz)
  {
    m_minx = minx;
    m_maxx = maxx;
    m_miny = miny;
    m_maxy = maxy;
    m_minz = minz;
    m_maxz = maxz;
  }
  float GetMinX() { return m_minx; }
  float GetMaxX() { return m_maxx; }
  float GetMinY() { return m_miny; }
  float GetMaxY() { return m_maxy; }
  float GetMinZ() { return m_minz; }
  float GetMaxZ() { return m_maxz; }
  private float m_minx;
  float m_maxx;
  float m_miny;
  float m_maxy;
  float m_minz;
  float m_maxz;
};

class Ball
{
  Ball(Coordinat coordinat, Vector speed, float radius, Arena arena)
  {
    m_button_state = 0.0;
    m_coordinat = coordinat;
    m_speed = speed;
    m_radius = radius;
    m_arena = arena;
  }
  void ChangeDy(float ddy) { SetDy(GetDy() + ddy); }
  void ChangeDz(float ddz) { SetDz(GetDz() + ddz); }
  void Display() 
  {
    translate(GetX(), GetY(), GetZ());
    sphere(GetRadius());
    translate(-GetX(), -GetY(), -GetZ());
  }
  void Move()
  {
    SetX(GetX() + GetDx());
    if (GetX() < GetMinX()) { SetDx( abs(GetDx())); }
    if (GetX() > GetMaxX()) { SetDx(-abs(GetDx())); }

    SetY(GetY() + GetDy());
    if (GetY() < GetMinY()) { SetDy( abs(GetDy())); }
    if (GetY() > GetMaxY()) { SetDy(-abs(GetDy())); }

    SetZ(GetZ() + GetDz());
    if (GetZ() < GetMinZ()) { SetDz( abs(GetDz())); }
    if (GetZ() > GetMaxZ()) { SetDz(-abs(GetDz())); }
  }
  float GetDx() { return m_speed.GetDx(); } 
  float GetDy() { return m_speed.GetDy(); } 
  float GetDz() { return m_speed.GetDz(); } 
  float GetX() { return m_coordinat.GetX(); } 
  float GetY() { return m_coordinat.GetY(); } 
  float GetZ() { return m_coordinat.GetZ(); } 
  float GetMinX() { return m_arena.GetMinX(); }
  float GetMaxX() { return m_arena.GetMaxX(); }
  float GetMinY() { return m_arena.GetMinY(); }
  float GetMaxY() { return m_arena.GetMaxY(); }
  float GetMinZ() { return m_arena.GetMinZ(); }
  float GetMaxZ() { return m_arena.GetMaxZ(); }
  float GetRadius() { return m_radius; }
  void SetButton(float state) { m_button_state = state; }
  void SetDx(float dx) { m_speed.SetDx(dx); }
  void SetDy(float dy) { m_speed.SetDy(dy); }
  void SetDz(float dz) { m_speed.SetDz(dz); }
  void SetX(float x) { m_coordinat.SetX(x); }
  void SetY(float y) { m_coordinat.SetY(y); }
  void SetZ(float z) { m_coordinat.SetZ(z); }
  Vector m_speed;
  Coordinat m_coordinat;
  float m_radius;
  Arena m_arena;
  float m_button_state;
}

//import processing.sound.*;
//SoundFile file;

final float speed_z = 0.01;
//final float speed_z = 0.001;

final float player1_r = 50;
final float player2_r = 50;
final float bal_r = 50.0;
final float minx = -500.0;
final float maxx =  500.0;
final float miny = -250.0;
final float maxy =  250.0;
final float minz = -250.0;
final float maxz =  250.0;

final float bal_x = 0.0;
final float bal_y = 0.0;
final float bal_z = 0.0;
final float bal_dx = 5.0;
final float bal_dy = 0.1;
final float bal_dz = 1.0;
final float player1_x = minx + (1.0 * player1_r);
final float player1_y = 0.0;
final float player1_z = 0.0;
final float player1_dy = 0.0;
final float player2_x = maxx - (1.0 * player2_r);
final float player2_y = 0.0;
final float player2_z = 0.0;
final float player2_dy = 0.0;
final float text_size = 128;
float hoek_z = 0.0;
int score1 = 0;
int score2 = 0;

Ball player_1 = new Ball(
  new Coordinat(player1_x, player1_y, player1_z), 
  new Vector(0.0, 0.0, 0.0), 
  player1_r,
  new Arena(minx, maxx, miny, maxy, minz, maxz)
);

Ball player_2 = new Ball(
  new Coordinat(player2_x, player2_y, player2_z), 
  new Vector(0.0, 0.0, 0.0), 
  player2_r,
  new Arena(minx, maxx, miny, maxy, minz, maxz)
);

Ball bal = new Ball(
  new Coordinat(bal_x, bal_y, bal_z), 
  new Vector(bal_dx, bal_dy, bal_dz), 
  bal_r,
  new Arena(minx, maxx, miny, maxy, minz, maxz) 
);

//Music
AudioPlayer player;
Minim minim;


Serial serial_port;



/// standalone means without Arduino
final boolean standalone = true;

void setup()
{
  fullScreen(P3D);
  textSize(text_size);
  if (!standalone)
  {
    serial_port = new Serial(this, Serial.list()[0], 9600);
  }

  minim = new Minim(this);
  player = minim.loadFile("EddiesTwister.mp3");
  player.play();
}

void draw()
{
  if (!standalone)
  {
    while (serial_port.available() > 0) 
    {
      final int spacer_value = serial_port.read();
      if (spacer_value != 0) continue; 
      if (serial_port.available() == 0) break;
      player_1.SetY(map(serial_port.read(), 1, 255, miny, maxy));
      if (serial_port.available() == 0) break;
      player_1.SetZ(map(serial_port.read(), 1, 255, minz, maxz));
      if (serial_port.available() == 0) break;
      player_1.SetButton(map(serial_port.read(), 1, 255, 0, 255));
      if (serial_port.available() == 0) break;
      player_2.SetY(map(serial_port.read(), 1, 255, miny, maxy));
      if (serial_port.available() == 0) break;
      player_2.SetZ(map(serial_port.read(), 1, 255, minz, maxz));
      if (serial_port.available() == 0) break;
      player_2.SetButton(map(serial_port.read(), 1, 255, 0, 255));
    }
  }
  background(32);
  ambientLight(64, 64, 64);
  noStroke();
  directionalLight(255, 0, 0, -cos(hoek_z * 5.0), -sin(hoek_z * 5.0), 0);
  directionalLight(0, 255, 0, 0, -cos(hoek_z * 7.0), -sin(hoek_z * 7.0));
  directionalLight(0, 0, 255, -sin(hoek_z * 3.0), 0, -cos(hoek_z * 3.0));

  translate(width / 2, height / 2);
  rotateX(hoek_z * 0.2);
  rotateY(hoek_z * 0.3);
  rotateZ(hoek_z);

  player_1.Move();
  player_2.Move();
  bal.Move();

  //bal.check_collision(player_1);
  //bal.check_collision(player_2);

  player_1.Display();
  player_2.Display();
  bal.Display();

  DrawArena();

  fill(255, 128, 128);
  text(score1, -text_size, text_size / 2);

  fill(128, 128, 255);
  text(score2,  text_size, text_size / 2);

  hoek_z += speed_z;
}

void DrawBlockLine(Coordinat from, Coordinat to)
{
  for (float i = 0; i < 100.0; i += 1.0)
  {
    final float block_x = from.GetX() + ((i / 100.0) * (to.GetX() - from.GetX()));
    final float block_y = from.GetY() + ((i / 100.0) * (to.GetY() - from.GetY()));
    final float block_z = from.GetZ() + ((i / 100.0) * (to.GetZ() - from.GetZ()));
    translate(block_x, block_y, block_z);
    box(10);
    translate(-block_x, -block_y, -block_z);
  }
}

void DrawArena()
{
  //Lines in length of the arena first
  DrawBlockLine(new Coordinat(minx, miny, minz), new Coordinat(maxx, miny, minz));
  DrawBlockLine(new Coordinat(minx, miny, maxz), new Coordinat(maxx, miny, maxz));
  DrawBlockLine(new Coordinat(minx, maxy, minz), new Coordinat(maxx, maxy, minz));
  DrawBlockLine(new Coordinat(minx, maxy, maxz), new Coordinat(maxx, maxy, maxz));

  //Lines at minx, just go head-tail
  DrawBlockLine(new Coordinat(minx, miny, minz), new Coordinat(minx, miny, maxz));
  DrawBlockLine(new Coordinat(minx, miny, maxz), new Coordinat(minx, maxy, maxz));
  DrawBlockLine(new Coordinat(minx, maxy, maxz), new Coordinat(minx, maxy, minz));
  DrawBlockLine(new Coordinat(minx, maxy, minz), new Coordinat(minx, miny, minz));

  //Lines at maxx, just go head-tail
  DrawBlockLine(new Coordinat(maxx, miny, minz), new Coordinat(maxx, miny, maxz));
  DrawBlockLine(new Coordinat(maxx, miny, maxz), new Coordinat(maxx, maxy, maxz));
  DrawBlockLine(new Coordinat(maxx, maxy, maxz), new Coordinat(maxx, maxy, minz));
  DrawBlockLine(new Coordinat(maxx, maxy, minz), new Coordinat(maxx, miny, minz));
}

float get_distance(float dx, float dy)
{
  return sqrt((dx * dx) + (dy * dy));
}

void keyPressed() {

  if (key == CODED) 
  {
    if (keyCode == UP) 
    {
      player_2.ChangeDy(-1);
    } 
    if (keyCode == DOWN) 
    {
      player_2.ChangeDy(1);
    }
    if (keyCode == LEFT) 
    {
      player_2.ChangeDz(-1);
    } 
    if (keyCode == RIGHT) 
    {
      player_2.ChangeDz(1);
    }
  }
  if (key == 'w') 
  {
    player_1.ChangeDy(-1);
  }
  if (key == 's') 
  {
    player_1.ChangeDy( 1);
  }
  if (key == 'a') 
  {
    player_1.ChangeDz(-1);
  }
  if (key == 'd') 
  {
    player_1.ChangeDz( 1);
  }
}

// From http://www.plasmaphysics.org.uk/programs/coll3d_cpp.htm:
/*****************************************************************************
//   This program is a 'remote' 3D-collision detector for two balls on linear
//   trajectories and returns, if applicable, the location of the collision for 
//   both balls as well as the new velocity vectors (assuming a partially elastic
//   collision as defined by the restitution coefficient).
//
//   All variables apart from 'error' are of Double Precision Floating Point type.
//
//   The Parameters are:
//
//    R    (restitution coefficient)  between 0 and 1 (1=perfectly elastic collision)
//    m1    (mass of ball 1)
//    m2    (mass of ball 2)
//    r1    (radius of ball 1)
//    r2    (radius of ball 2)
//  & x1    (x-coordinate of ball 1) 
//  & y1    (y-coordinate of ball 1)          
//  & z1    (z-coordinate of ball 1) 
//  & x2    (x-coordinate of ball 2)              
//  & y2    (y-coordinate of ball 2)         
//  & z2    (z-coordinate of ball 2)         
//  & vx1   (velocity x-component of ball 1) 
//  & vy1   (velocity y-component of ball 1)
//  & vz1   (velocity z-component of ball 1)          
//  & vx2   (velocity x-component of ball 2)         
//  & vy2   (velocity y-component of ball 2)
//  & vz2   (velocity z-component of ball 2)
//  & error (int)     (0: no error 
//                     1: balls do not collide
//                     2: initial positions impossible (balls overlap))
//
//   Note that the parameters with an ampersand (&) are passed by reference,
//   i.e. the corresponding arguments in the calling program will be updated 
//   (the positions and velocities however only if 'error'=0).
//   All variables should have the same data types in the calling program
//   and all should be initialized before calling the function.
//
//   This program is free to use for everybody. However, you use it at your own
//   risk and I do not accept any liability resulting from incorrect behaviour.
//   I have tested the program for numerous cases and I could not see anything 
//   wrong with it but I can not guarantee that it is bug-free under any 
//   circumstances.
//
//   I would appreciate if you could report any problems to me
//   (for contact details see  http://www.plasmaphysics.org.uk/feedback.htm ).
//
//   Thomas Smid   February 2004
//                 December 2005 (a few minor changes to improve speed)
//                 December 2009 (generalization to partially inelastic collisions)
//                 July     2011 (fix for possible rounding errors)
//******************************************************************************

   
    void collision3D(double R, double m1, double m2, double r1, double r2,
                     double& x1, double& y1,double& z1,
                     double& x2, double& y2, double& z2,
                     double& vx1, double& vy1, double& vz1,
                     double& vx2, double& vy2, double& vz2,
                     int& error)     {


       double  pi,r12,m21,d,v,theta2,phi2,st,ct,sp,cp,vx1r,vy1r,vz1r,fvz1r,
             thetav,phiv,dr,alpha,beta,sbeta,cbeta,dc,sqs,t,a,dvz2,
         vx2r,vy2r,vz2r,x21,y21,z21,vx21,vy21,vz21,vx_cm,vy_cm,vz_cm;

//     **** initialize some variables ****
       pi=acos(-1.0E0);
       error=0;
       r12=r1+r2;
       m21=m2/m1;
       x21=x2-x1;
       y21=y2-y1;
       z21=z2-z1;
       vx21=vx2-vx1;
       vy21=vy2-vy1;
       vz21=vz2-vz1;
       
       vx_cm = (m1*vx1+m2*vx2)/(m1+m2) ;
       vy_cm = (m1*vy1+m2*vy2)/(m1+m2) ;
       vz_cm = (m1*vz1+m2*vz2)/(m1+m2) ;  

     
//     **** calculate relative distance and relative speed ***
       d=sqrt(x21*x21 +y21*y21 +z21*z21);
       v=sqrt(vx21*vx21 +vy21*vy21 +vz21*vz21);
       
//     **** return if distance between balls smaller than sum of radii ****
       if (d<r12) {error=2; return;}
       
//     **** return if relative speed = 0 ****
       if (v==0) {error=1; return;}
       

//     **** shift coordinate system so that ball 1 is at the origin ***
       x2=x21;
       y2=y21;
       z2=z21;
       
//     **** boost coordinate system so that ball 2 is resting ***
       vx1=-vx21;
       vy1=-vy21;
       vz1=-vz21;

//     **** find the polar coordinates of the location of ball 2 ***
       theta2=acos(z2/d);
       if (x2==0 && y2==0) {phi2=0;} else {phi2=atan2(y2,x2);}
       st=sin(theta2);
       ct=cos(theta2);
       sp=sin(phi2);
       cp=cos(phi2);


//     **** express the velocity vector of ball 1 in a rotated coordinate
//          system where ball 2 lies on the z-axis ******
       vx1r=ct*cp*vx1+ct*sp*vy1-st*vz1;
       vy1r=cp*vy1-sp*vx1;
       vz1r=st*cp*vx1+st*sp*vy1+ct*vz1;
       fvz1r = vz1r/v ;
       if (fvz1r>1) {fvz1r=1;}   // fix for possible rounding errors
          else if (fvz1r<-1) {fvz1r=-1;} 
       thetav=acos(fvz1r);
       if (vx1r==0 && vy1r==0) {phiv=0;} else {phiv=atan2(vy1r,vx1r);}

                    
//     **** calculate the normalized impact parameter ***
       dr=d*sin(thetav)/r12;


//     **** return old positions and velocities if balls do not collide ***
       if (thetav>pi/2 || fabs(dr)>1) {
           x2=x2+x1;
           y2=y2+y1;
           z2=z2+z1;
           vx1=vx1+vx2;
           vy1=vy1+vy2;
           vz1=vz1+vz2;
           error=1;
           return;
        }
       
//     **** calculate impact angles if balls do collide ***
       alpha=asin(-dr);
       beta=phiv;
       sbeta=sin(beta);
       cbeta=cos(beta);
        
       
//     **** calculate time to collision ***
       t=(d*cos(thetav) -r12*sqrt(1-dr*dr))/v;

     
//     **** update positions and reverse the coordinate shift ***
       x2=x2+vx2*t +x1;
       y2=y2+vy2*t +y1;
       z2=z2+vz2*t +z1;
       x1=(vx1+vx2)*t +x1;
       y1=(vy1+vy2)*t +y1;
       z1=(vz1+vz2)*t +z1;
        
 
       
//  ***  update velocities ***

       a=tan(thetav+alpha);

       dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));
       
       vz2r=dvz2;
       vx2r=a*cbeta*dvz2;
       vy2r=a*sbeta*dvz2;
       vz1r=vz1r-m21*vz2r;
       vx1r=vx1r-m21*vx2r;
       vy1r=vy1r-m21*vy2r;

       
//     **** rotate the velocity vectors back and add the initial velocity
//           vector of ball 2 to retrieve the original coordinate system ****
                     
       vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
       vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
       vz1=ct*vz1r-st*vx1r               +vz2;
       vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
       vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
       vz2=ct*vz2r-st*vx2r               +vz2;
        

//     ***  velocity correction for inelastic collisions ***

       vx1=(vx1-vx_cm)*R + vx_cm;
       vy1=(vy1-vy_cm)*R + vy_cm;
       vz1=(vz1-vz_cm)*R + vz_cm;
       vx2=(vx2-vx_cm)*R + vx_cm;
       vy2=(vy2-vy_cm)*R + vy_cm;
       vz2=(vz2-vz_cm)*R + vz_cm;  

       return;
}

*/