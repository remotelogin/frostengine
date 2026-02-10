#define FILL_COLOR 0xAA22BB
#define BORDER_COLOR 0x000000
#define GLOBAL_ILLUMINATION_VALUE 0.1f

#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/XShm.h>
#include <iostream>
#include <vector>
#include "math.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <list>
#include <sys/shm.h>
#include <sys/ipc.h>

#include "utility.h"
#include "model.h"
#include "light.h"

void multiply_matrix_vector(vec3d &i, vec3d &o, mat4x4 &m) {

  o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
  o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
  o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];

  float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
  
  if (w != 0.0f) {

    o.x /= w; o.y /= w; o.z /= w;

  }

  return;
}

class XlibApp {
  
public:
  
  int sgn(int x)
  {
    if (x > 0) return +1;
    if (x < 0) return -1;
    return 0;
  }

  //this is a nightmare
  void draw_bresenham_line(unsigned int start_x, unsigned int start_y, unsigned int end_x, unsigned int end_y, unsigned long color) {

    XSetForeground(display,bb_gc,color);
    
    int x, y, t, dx, dy, incx, incy, pdx, pdy, ddx, ddy, deltaslowdirection, deltafastdirection, err;
    
    dx = end_x - start_x;
    dy = end_y - start_y;
    
    incx = sgn(dx);
    incy = sgn(dy);

    if (dx < 0) dx = -dx;
    if (dy < 0) dy = -dy;
    
    if (dx > dy) {
      
        pdx = incx; pdy = 0;
        ddx = incx; ddy = incy;
        deltaslowdirection = dy;
	deltafastdirection = dx;

    } else {

      pdx = 0;    pdy = incy;
      ddx = incx; ddy = incy;
      deltaslowdirection = dx;
      deltafastdirection = dy;   

    }
    
    x = start_x;
    y = start_y;
    err = deltafastdirection / 2;

    XDrawPoint(display,backBuffer,bb_gc,x,y);
    
    for(t = 0; t < deltafastdirection; ++t)  {

      err -= deltaslowdirection;

      if(err < 0) {

	err += deltafastdirection;
	x += ddx;
	y += ddy;

      } else {

	x += pdx;
	y += pdy;

      }

      XDrawPoint(display,backBuffer,bb_gc,x,y);

    }
  }

  // i mog nimma :(
  void draw_filled_in_triangle (unsigned int x1, unsigned int x2, unsigned int x3, unsigned int y1, unsigned int y2, unsigned int y3 ,unsigned long color) {

    float fx1 = x1;
    float fx2 = x2;
    float fx3 = x3;
    float fy1 = y1;
    float fy2 = y2;
    float fy3 = y3;
    
    std::cout << "drawing filled tri" << std::endl;
    
    unsigned int s_buf;
    
    if(fy1 > fy2) { s_buf = fy1; fy1 = fy2 ; fy2 = s_buf; }
    if(fy2 > fy3) { s_buf = fy2; fy2 = fy3 ; fy3 = s_buf; }
    
    if (fy2 == fy1)
      {
	std::cout << "bottom triangle! " << std::endl;
	
	fill_bottom_triangle(fx1,fy1,fx2,fy2,fx3,fy3,color);
      }
    
    if (fy2 == fy3)
      {
	std::cout << "top triangle! " << std::endl;
		
	fill_top_triangle(fx1,fy1,fx2,fy2,fx3,fy3,color);
      }
    
    float fx4,fy4;

    //    if(x1 == 0.0f || x2 == 0.0f || x3 == 0.0f || y1 == 0.0f || y2 == 0.0f || y3 == 0.0f) 
    
    fx4 = (fx1 + ((float)(fy2 - fy1) / (float)(fy3 - fy1)) * (fx3 - fx1));
    fy4 = fy2;

    fill_bottom_triangle(fx1,fy1,fx4,fy4,fx3,fy3,color);
    fill_top_triangle(fx1,fy1,fx4,fy4,fx3,fy3,color);
    
  }

  void draw_horizontal_line(unsigned int x1, unsigned int x2, unsigned int y1, unsigned int y2, unsigned long color) {

    XSetForeground(display,bb_gc,color);
    
    XDrawLine(display,backBuffer,bb_gc,x1,y1,x2,y2);

  }

  mesh import_obj_mesh(std::string file_path) {

    mesh output;
    
    std::ifstream file(file_path);
    if(!file.is_open())
      std::cout << "imort mesh not found!!!" << std::endl;

    std::vector<vec3d> vertices;
    
    while(!file.eof()) {

      char line[1024];
      file.getline(line,1024);

      std::stringstream s;

      s << line;

      char junk;

      if(line[0] == 'v') {

	vec3d vec;

	s >> junk >> vec.x >> vec.y >> vec.z;

        vertices.push_back(vec);

      }
      
      if(line[0] == 'f') {

	int f[3];

	s >> junk >> f[0] >> f[1] >> f[2];
	
	output.tris.push_back({ vertices[f[0] - 1], vertices[f[1] - 1], vertices[f[2] - 1] , 0.0f, 1.0f, 0x000000 });

      }
      
    }
    
    return output;

  }

  void try_to_render_screen() {    

    try_to_draw = true;

  }

  void window_runtime_helper() {

    printf("runtime helper start! \n");

    auto start_time = std::chrono::steady_clock::now();
    
    while(true) {            
      
      auto cur_time = std::chrono::steady_clock::now();

      auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(cur_time - start_time).count();

      avg_fps = (float)frames_drawn / (float)elapsed_seconds;
      
      if (try_move_foreward) {
	
	vec3d cam_foreward = vector_Mul(look_dir, 0.04f * z_pos);
	
	camera = vector_Add(camera, cam_foreward);	
      }
      
      if (try_move_backward) {
	
	vec3d cam_foreward = vector_Mul(look_dir, 0.04f * z_pos);
	
	camera = vector_Sub(camera, cam_foreward);	
      }

      if (try_move_right) {
	
	mat4x4 rot_x = matrix_make_rotate_y(-1.5707963268f);

	vec3d normal_look_dir = vector_Mul(look_dir,0.04f);

	vector_Normalise(normal_look_dir);
	
	vec3d cam_modified = matrix_multiply_vector(rot_x,normal_look_dir);
		
	camera = vector_Sub(camera, cam_modified);
	
	x_pos -= 0.01f;	 	
      }

      if (try_move_left) {
	
	mat4x4 rot_x = matrix_make_rotate_y(1.5707963268f);
	
	vec3d normal_look_dir = vector_Mul(look_dir,0.04f);

	vector_Normalise(normal_look_dir);
	
	vec3d cam_modified = matrix_multiply_vector(rot_x,normal_look_dir);
	
	camera = vector_Sub(camera, cam_modified);
	
	x_pos += 0.01f;	 	  	
      }
      
      if (try_move_down) 	
	camera.y -= 0.05f;	 
	      
      if (try_move_up) 	
	camera.y += 0.05f;	 
	     
      if (try_rotate_right) 	
	cam_yaw += change_rate * 0.2f;	 
	
      if (try_rotate_left) 	
	cam_yaw -= change_rate * 0.2f;	 
        
      if(draw_cooldown) {
	std::this_thread::sleep_for(std::chrono::milliseconds(50));
	draw_cooldown = false;
	
      }

      	render_screen();

	//+++++++++++++++++++++++++++

	loaded_lights[0].fPosition_X = std::sin(time/100) * 4;
	loaded_models[2].fPosition_X = loaded_lights[0].fPosition_X;
	
	loaded_lights[1].fPosition_X = std::sin(time2/80) * 3;
	loaded_models[3].fPosition_X = loaded_lights[1].fPosition_X;

	time2++;
	
	time++;
	
	//+++++++++++++++++++++++++++
	
      if(!draw_cooldown && try_to_draw) {
	
        draw_cooldown = true;
	try_to_draw = false;
	
	render_screen();

	draw_cooldown = false;
	
      } else {

	
	//std::this_thread::sleep_for(std::chrono::milliseconds(5));
	//	frames_drawn++;
	
      }
      
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
      
      if((try_move_backward || try_move_down || try_move_up || try_move_foreward || try_move_left || try_move_right || try_rotate_left || try_rotate_right))
	{
	  render_screen();
	}

    }
    printf("window runtime helper shutdown!\n"); 
  }
  
  void tmp_draw_filled_tri(unsigned int v1x, unsigned int v2x, unsigned int v3x, unsigned int v1y, unsigned int v2y, unsigned int v3y, unsigned long color) {
    
    XSetForeground(display,bb_gc,color);
    
    XPoint triangle[3];
    triangle[0].x = v1x; triangle[0].y = v1y;
    triangle[1].x = v2x; triangle[1].y = v2y;
    triangle[2].x = v3x; triangle[2].y = v3y;
    
    XFillPolygon(display, backBuffer ,bb_gc, triangle, 3, 2, 0);

  }
  
  void fill_top_triangle(unsigned int v1x, unsigned int v1y, unsigned int v2x, unsigned int v2y, unsigned int v3x, unsigned int v3y,unsigned int color) {

    float fv1x = (float) v1x;
    float fv1y = (float) v1y;
    float fv2x = (float) v2x;
    float fv2y = (float) v2y;    
    float fv3x = (float) v3x;
    float fv3y = (float) v3y;
    
    float invslope1,invslope2;
    
    if (fv3y != fv1y) {
      invslope1 = (fv3x - fv1x) / (fv3y - fv1y);
    } else {
      fv1y = fv1y + 0.1f;
      invslope1 = (fv3x - fv1x) / (fv3y - fv1y);
    }
    
    if (fv3y != fv2y) {
      invslope2 = (fv3x - fv2x) / (fv3y - fv2y);
    } else {
      fv1y = fv1y + 0.1f;
      invslope2 = (fv3x - fv2x) / (fv3y - fv2y);
    }
        
    float curx1 = fv3x;
    float curx2 = fv3x;
    
    for (int scanlineY = fv3y; scanlineY > fv1y; scanlineY--)
      {
	
	draw_horizontal_line((unsigned int)curx1, (unsigned int)scanlineY, (unsigned int)curx2, (unsigned int)scanlineY, color);
	curx1 -= invslope1;
	curx2 -= invslope2;
	
      }
  }
  
  void fill_bottom_triangle(unsigned int v1x,unsigned int v1y,unsigned int v2x,unsigned int v2y,unsigned int v3x,unsigned int v3y,unsigned int color) {

    float fv1x = (float) v1x;
    float fv1y = (float) v1y;
    float fv2x = (float) v2x;
    float fv2y = (float) v2y;    
    float fv3x = (float) v3x;
    float fv3y = (float) v3y;

    float invslope1,invslope2;
    if(fv2y != fv1y) {
      invslope1 = (fv2x - fv1x) / (fv2y - fv1y);
    } else {
      fv1y = fv1y + 0.1f;
      invslope1 = (fv2x - fv1x) / (fv2y - fv1y);
    }

    if(fv3y !=fv1y) {
      invslope2 = (fv3x - fv1x) / (fv3y - fv1y);
    } else {
      fv1y = fv1y + 0.1f;
      invslope2 = (fv3x - fv1x) / (fv3y - fv1y);
    }
    
    float curx1 = fv1x;
    float curx2 = fv1x;
    
    for (int scanlineY = fv1y; scanlineY <= fv2y; scanlineY++)
      {	
	draw_horizontal_line((unsigned int)curx1, (unsigned int)scanlineY, (unsigned int)curx2, (unsigned int)scanlineY, color);
	curx1 += invslope1;
	curx2 += invslope2;
      }
  }
  
  XlibApp(int width, int height) : width(width), height(height) {
    display = XOpenDisplay(nullptr);
    if (display == nullptr) {
      std::cerr << "Cannot open display" << std::endl;
      exit(1);
    }
    
    root = DefaultRootWindow(display);
    
    screen = DefaultScreen(display);

    window = XCreateSimpleWindow(display,
				 root,
				 10,
				 10,
				 width,
				 height,
				 1,				 
				 BlackPixel(display, screen),
				 WhitePixel(display, screen)
				 );
    
    XSelectInput(display, window, ExposureMask | KeyPressMask | KeyReleaseMask);
    XMapWindow(display, window);
    
    gc = XCreateGC(display, window, 0, nullptr);

    backBuffer = XCreatePixmap(display, window, width, height, DefaultDepth(display, screen));
    
    bb_gc = XCreateGC(display,backBuffer,0,nullptr);
    
  }
  
  ~XlibApp() {
    //painge
    shutdown = true;
    shmdt(shminfo.shmaddr);
    shmctl(shminfo.shmid, IPC_RMID, 0);   
    XFreeGC(display, gc);
    XDestroyWindow(display, window);
    XCloseDisplay(display);
  }

  
  unsigned long construct_rgb_value(unsigned long r, unsigned long g, unsigned long b) {
    return (
	    static_cast<unsigned int>(r) << 16) | 
      (static_cast<unsigned int>(g) << 8) | 
      (static_cast<unsigned int>(b));
  }
  
  unsigned long float_to_rgb_grayscale(float value) {

    if (value < GLOBAL_ILLUMINATION_VALUE) value = GLOBAL_ILLUMINATION_VALUE;
    if (value > 1.0f) value = 1.0f;
       
    unsigned char grayValue = static_cast<unsigned char>(value * 255);

    unsigned long rgb = (255 << 24) | (grayValue << 16) | (grayValue << 8) | grayValue;

    return rgb;

  }
  
  void draw_triangle(int x1, int y1, int x2, int y2, int x3, int y3) {
    
    XSetForeground(display, bb_gc, BlackPixel(display, DefaultScreen(display)));
    
    draw_bresenham_line(x1, y1, x2, y2, BORDER_COLOR);
    draw_bresenham_line(x2, y2, x3, y3, BORDER_COLOR); 
    draw_bresenham_line(x3, y3, x1, y1, BORDER_COLOR);
    //    XFlush(display);
  }
  
  void draw_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    
    Colormap colormap = DefaultColormap(display, screen);
    XColor color;
    color.red = (r << 8);
    color.green = (g << 8);
    color.blue = (b << 8);
    XAllocColor(display, colormap, &color);
    
    XSetForeground(display, bb_gc, color.pixel);
    XDrawPoint(display, backBuffer, bb_gc, x, y);
    //XFlush(display);
    
  }
  
  void render_screen() {
        
    auto filter_time = std::chrono::high_resolution_clock::now();
    
    std::vector<triangle> v_want_to_draw;
    
    for(auto mesh_to_render : loaded_models) {
                  
      matRotZ = matrix_make_rotate_z(mesh_to_render.fThetaZ*(3.1415f/180.0f));
      matRotX = matrix_make_rotate_x(mesh_to_render.fThetaX*(3.1415f/180.0f));
      matRotY = matrix_make_rotate_y(mesh_to_render.fThetaY*(3.1415f/180.0f));
      
      matProj = matrix_make_projection(fFov, (float)height / (float)width, 0.1f, 100.0f);
      
      //tjmat
      
      matTrans = matrix_make_translate(mesh_to_render.fPosition_X,mesh_to_render.fPosition_Y,mesh_to_render.fPosition_Z);
      
      //matrix allocation
      matWorld = matrix_make_static_identity();
      matWorld = matrix_mul_matrix(matWorld,matRotX);
      matWorld = matrix_mul_matrix(matWorld,matRotY);
      matWorld = matrix_mul_matrix(matWorld,matRotZ);
      matWorld = matrix_mul_matrix(matWorld,matTrans);
      
      vec3d up = {0,1,0};
      vec3d target = {0,0,1};
      mat4x4 cam_rot = matrix_make_rotate_y(cam_yaw);
      look_dir = matrix_multiply_vector(cam_rot, target);
      target = vector_Add(camera, look_dir);
      mat4x4 camera_mat = matrix_PointAt(camera, target, up);
      
      mat4x4 view_mat = matrix_QuickInverse(camera_mat);
      
      for(auto tri : mesh_to_render.model_data.tris) {
	
	triangle triProjected; 
	triangle triTransformed;
	triangle triViewed;
	
	triTransformed.p[0] = matrix_multiply_vector(matWorld,tri.p[0]);
	triTransformed.p[1] = matrix_multiply_vector(matWorld,tri.p[1]);
	triTransformed.p[2] = matrix_multiply_vector(matWorld,tri.p[2]);
	
	vec3d normal, vec1, vec2;
	
	//vecs that span surface
	vec1 = vector_Sub(triTransformed.p[1], triTransformed.p[0]);
	vec2 = vector_Sub(triTransformed.p[2], triTransformed.p[0]);
	
	// normal of surface
	normal = vector_CrossProduct(vec1,vec2);
	normal = vector_Normalise(normal);
	
	vec3d camera_ray = vector_Sub(triTransformed.p[0], camera);
	
	if( vector_DotProduct(normal, camera_ray) < 0.0f ) {

	  /* not sure if this is shit impl but it is what it is
	  vec3d light_source = { 0.0f, 5.0f, -1.0f };
	  light_source = vector_Normalise(light_source);
	  float dp_light = normal.x * light_source.x + normal.y * light_source.y + normal.z * light_source.z;
	  */

	  //tjlight

	  float col_rgb_towrite = 0;
	  
	  for(auto i_light : loaded_lights) {

	    if(i_light.type == FACE_SHADE) { 

	      vec3d light_location = i_light.get_light_source();
	      
	      vector_Normalise(light_location);
	      
	      float delta_to_light = normal.x * i_light.fPosition_X + normal.y * i_light.fPosition_Y + normal.z * i_light.fPosition_Z;
	      
	      col_rgb_towrite = delta_to_light;
	      
	    }

	    if(i_light.type == FACE_SHADE_DISTANCE) {
	      
	      vec3d P = triTransformed.p[0];
	      float dx = P.x - i_light.fPosition_X;
	      float dy = P.y - i_light.fPosition_Y;
	      float dz = P.z - i_light.fPosition_Z;
	      float distance_to_light = std::sqrt(dx*dx + dy*dy + dz*dz);
	      
	      if(distance_to_light > i_light.strength){
		triTransformed.col_rgb  = float_to_rgb_grayscale( 0x000000 ); 
                continue;
	      }

	      float normalized_light_strength = 1-(0.001f + (( 0.999f - 0.001f ) / i_light.strength ) * distance_to_light); 

	      vec3d L = { i_light.fPosition_X - P.x,
			  i_light.fPosition_Y - P.y,
			  i_light.fPosition_Z - P.z };
	      vector_Normalise(L);
	      float ndotl = std::max(0.0f, vector_DotProduct(normal, L));
             
	      normalized_light_strength = normalized_light_strength * ndotl;
	      
	      col_rgb_towrite = col_rgb_towrite + normalized_light_strength;
	      
	    }

	  }

	  if (col_rgb_towrite > 0.99f)
	    col_rgb_towrite = 0.99f;
	  
	  triTransformed.col_rgb  = float_to_rgb_grayscale( col_rgb_towrite ); 
	  
	  // world -> view
	  triViewed.p[0] = matrix_multiply_vector(view_mat, triTransformed.p[0]);
	  triViewed.p[1] = matrix_multiply_vector(view_mat, triTransformed.p[1]);
	  triViewed.p[2] = matrix_multiply_vector(view_mat, triTransformed.p[2]);
	  
	  int nClippedTriangles = 0;
	  triangle clipped[2];
	  nClippedTriangles = clip_triangle_plane({ 0.0f, 0.0f, 0.5f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);
	  
	  for (int n = 0; n < nClippedTriangles; n++)
	    {
	      
	      // proj to 2d
	      triProjected.p[0] = matrix_multiply_vector(matProj, clipped[n].p[0]);
	      triProjected.p[1] = matrix_multiply_vector(matProj, clipped[n].p[1]);
	      triProjected.p[2] = matrix_multiply_vector(matProj, clipped[n].p[2]);
	      
	      triProjected.col_rgb = triTransformed.col_rgb;
	      // TJtrans
	      triProjected.p[0] = vector_Div(triProjected.p[0], triProjected.p[0].w);
	      triProjected.p[1] = vector_Div(triProjected.p[1], triProjected.p[1].w);
	      triProjected.p[2] = vector_Div(triProjected.p[2], triProjected.p[2].w);
	      
	      offset_view = {1,1,0};
	      triProjected.p[0] = vector_Add(triProjected.p[0], offset_view);
	      triProjected.p[1] = vector_Add(triProjected.p[1], offset_view);
	      triProjected.p[2] = vector_Add(triProjected.p[2], offset_view);
	      
	      triProjected.p[0].x *= 0.5f * (float)width;
	      triProjected.p[0].y *= 0.5f * (float)height;
	      
	      triProjected.p[1].x *= 0.5f * (float)width;
	      triProjected.p[1].y *= 0.5f * (float)height;
	      
	      triProjected.p[2].x *= 0.5f * (float)width;
	      triProjected.p[2].y *= 0.5f * (float)height;
	      
	      v_want_to_draw.push_back(triProjected);
	      
	    }
	  
	}      

      }
      
      //shoutout to stdlib
      sort(v_want_to_draw.begin(), v_want_to_draw.end(), [](triangle &t1, triangle &t2) {
	
	float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
	float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
	
	return z1 > z2;
	
      });      
      
    }

    auto render_start = std::chrono::high_resolution_clock::now();
    
    //done with each component    
    
    //clear canvas
    XSetForeground(display,bb_gc,0x000000);
    XFillRectangle(display, backBuffer, bb_gc, 0, 0, width, height);     
    
    for(auto &tri_rasterized : v_want_to_draw) {
      
      triangle clipped[2];
      std::list<triangle> listTriangles;
      listTriangles.push_back(tri_rasterized);
      int newTriangles = 1;

      for (int p = 0; p < 4; p++)
	{
	  int nTrisToAdd = 0;
	  while (newTriangles > 0)
	    {
	      
	      triangle test = listTriangles.front();
	      listTriangles.pop_front();
	      newTriangles--;
	      switch (p)
		{
		case 0:	nTrisToAdd = clip_triangle_plane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
		case 1:	nTrisToAdd = clip_triangle_plane({ 0.0f, (float)height - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
		case 2:	nTrisToAdd = clip_triangle_plane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
		case 3:	nTrisToAdd = clip_triangle_plane({ (float)width - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
		}
	      
	      for (int w = 0; w < nTrisToAdd; w++)
		listTriangles.push_back(clipped[w]);
	    }
	  newTriangles = listTriangles.size();
	}

      for(auto &t : listTriangles) {      
	
	tmp_draw_filled_tri((unsigned int)t.p[0].x,
			    (unsigned int)t.p[1].x,
			    (unsigned int)t.p[2].x,
			    (unsigned int)t.p[0].y,
			    (unsigned int)t.p[1].y,
			    (unsigned int)t.p[2].y,
			    t.col_rgb);
	
	if(wireframe_active){

	  //why r bresenhams so slow lol
          draw_triangle(t.p[0].x, t.p[0].y,
			t.p[1].x, t.p[1].y,
			t.p[2].x, t.p[2].y);
	  
	}      
	
      }
      
    }

    auto render_end = std::chrono::high_resolution_clock::now();
    
    if(wireframe_active) {
      XSetForeground(display, bb_gc, 0x00FF00);
      XDrawString(display, backBuffer , bb_gc, 10, 50, "Wireframe - active", 18);
    } else {
      XSetForeground(display, bb_gc, 0xFF0000);
      XDrawString(display, backBuffer , bb_gc, 10, 50, "Wireframe - deactivated", 23);
    }

    auto frame_time = std::chrono::duration_cast<std::chrono::nanoseconds>(render_end - render_start);
    auto filter_time_calc = std::chrono::duration_cast<std::chrono::nanoseconds>(render_start - filter_time);
    
    std::string cur_fov_str = "fFov - ";
    cur_fov_str += std::to_string(fFov);
    XSetForeground(display, bb_gc, 0x00FF00);
    XDrawString(display, backBuffer , bb_gc, 10, 80, cur_fov_str.c_str(), cur_fov_str.length());

    std::string cur_fps_str = "FPS (average) - ";
    cur_fps_str += std::to_string(avg_fps);
    XSetForeground(display, bb_gc, 0xAAAAFF);
    XDrawString(display, backBuffer , bb_gc, 10, 110, cur_fps_str.c_str(), cur_fps_str.length());

    std::string filter_str = "vertex filter time (uS) - ";
    filter_str += std::to_string(filter_time_calc.count() / 1000);
    XSetForeground(display, bb_gc, 0xAAAAFF);
    XDrawString(display, backBuffer , bb_gc, 10, 140, filter_str.c_str(), filter_str.length());
    
    std::string frame_time_str = "Frame time (uS)- ";
    frame_time_str += std::to_string(frame_time.count() / 1000);
    XSetForeground(display, bb_gc, 0xAAAAFF);
    XDrawString(display, backBuffer , bb_gc, 10, 170, frame_time_str.c_str(), frame_time_str.length());
    
    std::string num_polygons = "loaded Polygons - ";
    num_polygons +=  std::to_string(loaded_models[0].model_data.tris.size());
    XSetForeground(display, bb_gc, 0xFFAAFF);
    XDrawString(display, backBuffer , bb_gc, 10, 200, num_polygons.c_str(), num_polygons.length());

    // swap buffers??? hah i wish 

    XSync(display, False);
    
    XCopyArea(display, backBuffer, window, gc, 0, 0, width, height, 0, 0);

    XSync(display, False);
	
    // SWAPPING BUFFERS!!!??
    
    frames_drawn++;   
    
    draw_cooldown = false;
    
  }    

  void load_model(std::string to_load, float pos_x_in, float pos_y_in, float pos_z_in, float rot_x_in, float rot_y_in, float rot_z_in  ) {

    model imported_model = model(import_obj_mesh(to_load) , pos_x_in , pos_y_in , pos_z_in , rot_x_in , rot_y_in , rot_z_in);
    
    loaded_models.push_back(imported_model);

    return;
  }

  void load_light(unsigned int strength,
                float pos_x_in, float pos_y_in, float pos_z_in,
                float rot_x_in, float rot_y_in, float rot_z_in,
                available_light_types type)
{
  model imported_model = model(import_obj_mesh("models/light_source_mesh.obj"),
                               pos_x_in, pos_y_in, pos_z_in,
                               rot_x_in, rot_y_in, rot_z_in);

  light imported_light = light(strength,
                               rot_x_in, rot_y_in, rot_z_in,
                               pos_x_in, pos_y_in, pos_z_in,
                               type);

  loaded_models.push_back(imported_model);
  loaded_lights.push_back(imported_light);
}

  
  
  void run() {

    printf("launching draw helper! \n");
    
    std::thread t(&XlibApp::window_runtime_helper,this);
    
    camera.x = 0.0f;
    camera.y = 0.0f;
    camera.z = 0.0f;
  
    screen = DefaultScreen(display);

    shminfo.shmid = shmget(IPC_PRIVATE, width * height * 4, IPC_CREAT | 0777);
    shminfo.shmaddr = (char*)shmat(shminfo.shmid, 0, 0);
    shminfo.readOnly = False;

    if (shminfo.shmaddr == (void *)-1) {
      std::cout << "unable to allocate shm!" << std::endl;
      return;
    }
    
    XEvent event;
    while (true) {

      // pixmaps suck :(

      XNextEvent(display, &event);
      
      if (event.type == Expose) {
	
	int x, y; 
	unsigned int border_width, depth;
	
	if (XGetGeometry(display, window, &root, &x, &y, &width, &height, &border_width, &depth)) {
	  std::cout << "Window size: " << width << "x" << height << std::endl;
	} else {
	  std::cerr << "Failed to get window geometry" << std::endl;
	}
	
        try_to_render_screen();
	
      }
      
      if (event.type == KeyPress) {
	KeySym key = XLookupKeysym(&event.xkey, 0);

	// world rotation
	if (key == XK_n) {
	  
	  fThetaX += change_rate;     
	  
	} else if (key == XK_k) {
	  
	  fThetaY -= change_rate;	  
	  
	} else if (key == XK_m) {
	  
	  fThetaX -= change_rate;	  
	  
	} else if (key == XK_j) {
	  
	  fThetaY += change_rate;	  
	  
	} else if (key == XK_u) {
	  
	  fThetaZ += change_rate;	  
	  
	} else if (key == XK_i) {
	  
	  fThetaZ -= change_rate;	 
	  
	} else if (key == XK_f) {
	  
	  wireframe_active = !wireframe_active;
	  
	}
	
	// camera navigation input

	if (key == XK_w)
	  try_move_foreward = true;

	if (key == XK_s)
	  try_move_backward = true;

	if (key == XK_a)
	  try_move_right = true;

	if (key == XK_d)
	  try_move_left = true;

	if (key == XK_h)
	  try_rotate_right = true;

	if (key == XK_l)
	  try_rotate_left = true;

	if (key == XK_space)
	  try_move_down = true;

	if (key == XK_Shift_L)
	  try_move_up = true;

	//extras

	if (key == XK_t) {
	  
	  fFov += 5.0f;
	  //fFovRad = 1.0f / tanf( fFov * 0.5f / 180.0f * 3.14159f );
	  fFovRad += 0.1f;
	}
	
	if (key == XK_g) {
	  
	  fFov -= 5.0f;
	  //fFovRad = 1.0f / tanf( fFov * 0.5f / 180.0f * 3.14159f );
	  fFovRad -= 0.1f;
	}
	
	//quit
	if (key == XK_x) {
	  break; 
	}
	
        try_to_render_screen();
	
      }

      if(event.type == KeyRelease) {
	
	KeySym key = XLookupKeysym(&event.xkey, 0);

	if (key == XK_w)
	  try_move_foreward = false;
	if (key == XK_s)
	  try_move_backward = false;
	if (key == XK_a)
	  try_move_right = false;
	if (key == XK_d)
	  try_move_left = false;
	if (key == XK_h)
	  try_rotate_right = false;
	if (key == XK_l)
	  try_rotate_left = false;
	if (key == XK_space)
	  try_move_down = false;
	if (key == XK_Shift_L)
	  try_move_up = false;

      }				      
      
    }
    
  }
  
  
private:
  
  //buffer ovrflow speedrun
  float time = 0;
  float time2 = 0;
  
  float change_rate = 0.15f;
  
  bool wireframe_active = false;

  //X11
  Display* display;
  Window window;
  GC gc;
  Window root;
  int screen;
  XShmSegmentInfo shminfo;
  
  Pixmap backBuffer;
  GC bb_gc;
  
  unsigned int width;
  unsigned int height;

  //Renderer Data
  float fThetaX = 0.0f;
  float fThetaY = 0.0f;
  float fThetaZ = 0.0f;

  float fPosition_X;
  float fPosition_Y;
  float fPosition_Z;
  
  mesh loaded_mesh;

  float z_pos = 1.0f;
  float x_pos = 1.0f;
  
  float fNear = 0.1f;
  float fFar = 1000.0f;
  float fFov = 100.0f;

  vec3d camera;
  vec3d look_dir;
  vec3d offset_view;
  float cam_yaw;
  
  mat4x4 matProj;
  mat4x4 matRotZ;
  mat4x4 matRotX;
  mat4x4 matRotY;

  mat4x4 matTrans;
  mat4x4 matWorld;
  
  float fAspectRatio = (float)height / (float)width;
  float fFovRad = 1.72123f;

  int clip_margin = 50;
  
  //Optimization and debug
  bool try_to_draw;
  bool draw_cooldown = false;
  bool shutdown = false;

  unsigned int frames_drawn = 0;
  float avg_fps = 0;

  //Movement helpers (i know its shit impl but its x11 so theres no better way )

  bool try_move_foreward = false;
  bool try_move_backward = false;
  bool try_move_left = false;
  bool try_move_right = false;
  bool try_move_up = false;
  bool try_move_down = false;
  bool try_rotate_left = false;
  bool try_rotate_right = false;

  //data structure to hold abstracted models
  std::vector<model> loaded_models;
  std::vector<light> loaded_lights;
  
};

int main() {
  //init window
  XlibApp app(1920, 1080);

  //load models
  app.load_model("models/keyboard.obj", 0.0f , 0.0f , 0.0f , 0.0f , 0.0f , 0.0f );
  app.load_model("models/floor.obj", 0.0f , 3.0f , 0.0f , 180.0f , 0.0f , 0.0f );

  app.load_light(5,-3.0f,0.0f,-3.0f,0.0f,0.0f,0.0f,FACE_SHADE_DISTANCE);
  app.load_light(5,-3.0f,0.0f,0.0f,0.0f,0.0f,0.0f,FACE_SHADE_DISTANCE);
  
  //start engine
  app.run();
  return 0;
}
