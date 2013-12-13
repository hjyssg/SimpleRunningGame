//
//  RunningGameController.cpp
//  Running
//
//  Created by Junyang Huang on 11/6/13.
//  Copyright (c) 2013 . All rights reserved.
//


#include "HJY_GL.cpp"
#include "RunningGameEngine.cpp"
#include "tiny_obj_loader.h"

#include <vector>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cassert>
#include <pthread.h>

#include "Texture.h"


void glut_setup(void) ;
void gl_setup(void) ;
void init_game();
void load_files();
void my_init(int argc, char **argv);
void my_display(void) ;
void my_reshape(int w, int h) ;
void my_mouse(int x, int y);
void my_mouse2(int button, int state, int x, int y);
void my_keyboard_up( unsigned char key, int x, int y ) ;
void my_keyboard( unsigned char key, int x, int y ) ;
void my_special_keyboard(int key, int x, int y);
void timeElapse(int value);
int LoadGLTextures();
void play_music(char *filename);
void *play_music_thread_method(void *filename);
bool LoadTGA(Texture *, char *);


const int windows_height = 900;
const int windows_width = 600;

const float street_height = 20;


bool display_instruction_page = false;

//a flag to stop the game
bool stop = false;

//decide if to turn on audio
bool music_on = true;

bool look_back = false;

//when hit a tree, a hurrting effect
//then human will flickering
int hurting_effect_counter = 0;

const int number_frame = 32;
const int frame_time=  1000/number_frame;
const float time_step = (float)1/number_frame;
int millsecond_counter = 0;

const int start_screen_time = frame_time* number_frame *4; //for second
//a countdown before the game
int start_counter = start_screen_time;


typedef enum {WIN_THE_GAME, LOSE_THE_GAME, GAME_ON}GAME_CONDITION;
GAME_CONDITION game_condition = GAME_ON;


GLfloat colors [][3] = {
    {0.1, 0.2, 0.34},
    {0.96, 0.0, 0.0},  /* red     */
    {1.0, 0.0, 1.0},  /* magenta */
    {0.0, 0.0, 1.0},  /* blue    */
    {0.5, 0.5, 0.5},  /* 50%grey */
    {1.0, 1.0, 1.0} ,  /* white   */
};

GLfloat tree_colors[][3] =
{
    {0.654, 0.9, 0.61}, //green 1
    {0.34, 0.8, 0.31}, //green 3
    {0.2, 0.8, 0.1} //green 3
};

GLdouble coin_angle = 0;





const int GOAL_DISTANCE  = 1200;
const int runner_num = 4;

Texture * ground_texture;
Texture *sky_texture;
Texture *instruction_page;

std::vector<Human *>humanList;  //the redenering human

std::vector<tinyobj::shape_t> tree;
std::vector<tinyobj::shape_t> tree2;
std::vector<tinyobj::shape_t> tree3;


#pragma mark input file path
char tree_file_name[] = "resource/cartoontrees/tree1.obj";
char tree2_file_name[] = "resource/cartoontrees/tree2.obj";
char tree3_file_name[] = "resource/cartoontrees/tree3.obj";



char ground_texture_filename[] = "resource/texture/forest_pathway.tga";
char sky_file_name[] = "resource/texture/sky_107.tga";
char instruction_fn[] = "resource/texture/instruction_page.tga";

char start_music_fn[] = "resource/sound/herewego.WAV";
char hit_tree_music_fn[] = "resource/sound/pain.WAV";
char coin_music_fn[] = "resource/sound/coin.WAV";
char lose_music_fn[] = "resource/sound/gameover.WAV";
char win_music_fn[] = "resource/sound/yahoo.WAV";

char game_title[] = "Simple Runnging Game";

//the camera object
HJYCamera * my_cam = HJYCamera::getCamera();

//represent the internal data of the game
int road_width = 10;
int view_boundry = 130;
RunningWorld * engine;
Player *user;



int main(int argc, char** argv)
{
     init_game();
    
    srand ( (unsigned int)time(NULL) );
    //RunningWorld::game_engine_unit_test();
    setbuf(stdout, NULL);   /* for writing to stdout asap */
    glutInit(&argc, argv);
    
    srand((unsigned)time(0));
    glut_setup();
    gl_setup();
    
    load_files();
    
    for (int ii = 0; ii < runner_num; ii++)
    {
        Human *h = new Human();
        humanList.push_back(h);
    }
    
    glutMainLoop();
    return(0);
}

void load_files()
{
    LoadGLTextures();
    
    std::string err2 = tinyobj::LoadObj(tree, tree_file_name, NULL);
    std::string err3 = tinyobj::LoadObj(tree2, tree2_file_name, NULL);
    std::string err4 = tinyobj::LoadObj(tree3, tree3_file_name, NULL);
}

void init_game()
{
    engine = new RunningWorld( view_boundry, road_width, runner_num-1);
    user = &engine->human_player;
    
    start_counter = start_screen_time;
    millsecond_counter = 0;
    stop = false;
    
}

#pragma mark glut-related
void glut_setup ()
{
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    
    glutInitWindowSize(windows_height,windows_width);
    glutInitWindowPosition(20,20);
    glutCreateWindow(game_title);
    
    /* set up callback functions */
    glutDisplayFunc(my_display);
    glutReshapeFunc(my_reshape);
    
    glutMouseFunc(my_mouse2);
    glutPassiveMotionFunc(my_mouse);
    glutMotionFunc(my_mouse);
    glutKeyboardFunc(my_keyboard);
    glutSpecialFunc(my_special_keyboard);
    
    timeElapse(0);
    
    return;
}

//bind data with opengl texure
//bacuase mose texuture of this project parameters are the same
//doing so to reduce code length
void bind_openg_texture(Texture * text, GLint internalformat, GLenum format)
{
    glGenTextures(1, &text->texID);
    glBindTexture(GL_TEXTURE_2D, text->texID);
    
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, internalformat, text->width, text->height, 0, format, GL_UNSIGNED_BYTE, text->imageData);
    
    if (text->imageData)
    {
        free(text->imageData);					// Free The Texture Image Memory ( CHANGE )
    }
    
}

int LoadGLTextures()
{
    ground_texture = (Texture *)malloc(sizeof(Texture));
    sky_texture = (Texture *) malloc(sizeof(Texture));
    instruction_page = (Texture *) malloc(sizeof(Texture));
    
    if (LoadTGA(ground_texture,ground_texture_filename))
    {
        bind_openg_texture(ground_texture, GL_RGB, GL_RGBA);
        printf("load ground texute\n");
    }
    
    
    
    if (LoadTGA(sky_texture, sky_file_name))
    {
        bind_openg_texture(sky_texture, GL_RGB, GL_RGB);
        printf("load street texute\n");
    }
    
    return 0;
}


bool orthengineode = false;
void set_projecttion()
{
    glMatrixMode(GL_PROJECTION) ;
    glLoadIdentity() ;
    my_cam->perspective();
    glMatrixMode(GL_MODELVIEW) ;
}

void gl_setup(void)
{
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH); //Gouraud
    
    // enable depth handling (z-buffer)
    glEnable(GL_DEPTH_TEST);
    
    // enable auto normalize
    glEnable(GL_NORMALIZE);
    
    // define the background color
    glClearColor(1,1,1,1);
    
    
    
    my_cam->position.x = 0;
    my_cam->position.y = 2;
    my_cam->position.z = 2;
    
    my_cam->lookAt.x = 0;
    my_cam->lookAt.y = 0;
    my_cam->lookAt.z = 0;
    
    my_cam->up.x = 0;
    my_cam->up.y = 1;
    my_cam->up.z = 0;
    my_cam->up.normalize();
    
    my_cam->near = 1.5;
    my_cam->far = view_boundry*4;
    my_cam->zoom(2.2);
    
    
    set_projecttion();
    
    
    // toggle to smooth shading (instead of flat)
    glShadeModel(GL_SMOOTH);
    
    return ;
}



void my_reshape(int w, int h)
{
    glViewport(0,0,w,h) ;
    return ;
}

void my_special_keyboard(int key, int x, int y)
{
    switch (key)
    {
        case GLUT_KEY_DOWN:
            my_keyboard('s', 0, 0);
            break;
        case GLUT_KEY_UP:
            my_keyboard('w', 0, 0);
            break;
        case GLUT_KEY_LEFT:
            my_keyboard('a', 0, 0);
            break;
        case GLUT_KEY_RIGHT:
            my_keyboard('d', 0, 0);
            break;
            
    }
}

void my_mouse2(int button, int state,
               int x, int y)
{
    switch( button )
    {
           
            
        case GLUT_LEFT_BUTTON:
            if( state == GLUT_DOWN )
            {
                
                user->spdUp_by_using_coin();
            }
            
            if( state == GLUT_UP )
            {
            }
            break ;
    }
    
}

void my_mouse(int x, int y)
{
    static int previous_x = x;
    const int sld = 1;
    static int slow_conter = sld;
    
    if (slow_conter == 0) {
     
        //printf("x %i, y %i\n", x, y);
        
        //only response when the cursor is inside the windows
        if (x<0||y<0||x>windows_width||y>windows_height)
        {
            return;
        }
        
        if (x < previous_x ) {
            
            engine->move_player(user, -1);
        }else
        {
            engine->move_player(user, 1);
        }
        
        slow_conter = sld;
    }
    
    slow_conter--;
}

void my_keyboard( unsigned char key, int x, int y )
{
    switch( key )
    {
        case ' ':{
            user->spdUp_by_using_coin();
        }break;
            
        case '0':{
            engine->player_jump(user);
        }break;
            
        case 'p':
        case 'P':
        {
            engine->debug_print();
        }
            break;
            
        case 's':
        case 'S':
        {
            
        }
            break;
        case 'w':
        case 'W':
        {
            
        }break;
            
            
        case 'A':
        case 'a':
        {
            
            engine->move_player(user,-1);
        }break;
        case 'd':
        case 'D':
        {
            engine->move_player(user,1);
            
        }break;
            
            
        case '1':
        {
            music_on = !music_on;
        }
            break;
            
            
        case '2':
        {
            //printf("stop");
            stop = !stop;
            
        }break;
            
        case '3':
        {
            look_back = !look_back;
        }break;
            
            
        case '+':
        {
            my_cam->zoom(0.9);
            set_projecttion();
            glutPostRedisplay();
            
            printf("0.9 ");
            
        }break;
            
            
        case '-':
        case '_':
        {
            my_cam->zoom(1.1);
            set_projecttion();
            glutPostRedisplay();
            
            printf("1.1 ");
            
        }break;
            
            
        case 'e':
        case 'E':
        {
            engine->auto_mode = !engine->auto_mode;
            
            if (engine->auto_mode)
            {
                glutSetWindowTitle("auto mode on");
            }else
            {
                glutSetWindowTitle("auto mode off");
            }
        }
            break;
            
        case 'q':
        case 'Q':
        {
            exit(0) ;
            break ;
        }
        default: break;
    }
    
    return ;
}

#pragma mark calculation

inline bool isInUserViewRange(float z)
{
        //look forward
        if (my_cam->lookAt.z < 0 )
        {
            return  z < (user->z_pos + my_cam->far) && z > user->z_pos;
        }else
        {
            return z > (user->z_pos - view_boundry) && z < user->z_pos;
        }
    
//    if (z > user->z_pos - view_boundry&& z < user->z_pos + my_cam->far)
//    {
//    
//        
//    }else
//    {
//        return false;
//    }
}

inline int get_elapsed_time()
{
    return millsecond_counter/1000;
}

inline float get_opengl_z(PathObject *obj)
{
    return -obj->z_pos;
}

inline int get_opengl_x(PathObject *obj)
{
    return (obj->x_pos - road_width/2);
}
#pragma mark drawing method
void draw_str(char *str)
{
    
    
    for (int ii = 0; ii< strlen(str); ii++)
    {
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, str[ii]);
    }
}

void display_obj_mesh( std::vector<tinyobj::shape_t> shapes)
{
    
    //every three index form a face
    //their value reflect to vertex directly
    for (size_t i = 0; i < shapes.size()&&i <40; i++)
    {
#if 0
        tinyobj::material_t * mat = &shapes[i].material;
        
        if (mat->ambient[0]!=0.0 && mat->ambient[1]!=0.0 &&mat->ambient[2]!=0.0)
        {
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat->ambient);
        }
        
        if (mat->diffuse[0]!=0.0 && mat->diffuse[1]!=0.0 &&mat->diffuse[2]!=0.0)
        {
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat->diffuse);
        }
        
        if (mat->specular[0]!=0.0 && mat->specular[1]!=0.0 &&mat->specular[2]!=0.0)
        {
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat->shininess);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat->specular);
        }
#endif
        
        
        for (size_t f = 0; f < shapes[i].mesh.indices.size(); f += 3)
        {
            
            float v[3][3];
            float n[3][3];
            float t[3][2];
            
            for ( int  ii = 0; ii<3; ii++)
            {
                
                for (int jj = 0; jj<3; jj++)
                {
                    int id = shapes[i].mesh.indices[f+ii];
                    v[ii][jj] = shapes[i].mesh.positions[3*id+jj];
                    n[ii][jj] = shapes[i].mesh.normals[3*id+jj];
                }
            }
            
            if (shapes[i].mesh.texcoords.size()>0)
            {
                for (int ii = 0; ii<3; ii++)
                {
                    int id = shapes[i].mesh.indices[f+ii];
                    t[ii][0] = shapes[i].mesh.texcoords[2*id];
                    t[ii][1] = shapes[i].mesh.texcoords[2*id+1];
                }
            }
            
            glBegin(GL_POLYGON);
            glNormal3fv(n[0]);
            glVertex3fv(v[0]);
            glNormal3fv(n[1]);
            glVertex3fv(v[1]);
            glNormal3fv(n[2]);
            glVertex3fv(v[2]);
            glEnd();
        }
        
#if 0
        for (size_t v = 0; v < shapes[i].mesh.normals.size() / 3; v++) {
            printf("  n[%ld] = (%f, %f, %f)\n", v,
                   shapes[i].mesh.normals[3*v+0],
                   shapes[i].mesh.normals[3*v+1],
                   shapes[i].mesh.normals[3*v+2]);
        }
#endif
        
    }
}



void draw_tree(SlowDownTrap *tree_obj)
{
    glPushMatrix();
    {
        
        static float y_offset = -0.5;
        
        glTranslatef(get_opengl_x(tree_obj), y_offset, get_opengl_z(tree_obj));
        glScalef(0.17, 0.7, 0.17);
        glColor3f(0.1, 0.1, 0.1);
        GLfloat mat_shininess[] = { 50.0 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, tree_colors[tree_obj->type]);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tree_colors[tree_obj->type]);
        
        
        if (tree_obj->type == SLOWDOWN_TYPE_1)
        {
            display_obj_mesh(tree);
        }else if (tree_obj->type == SLOWDOWN_TYPE_2)
        {
            display_obj_mesh(tree2);
        }
        else if (tree_obj->type == SLOWDOWN_TYPE_3)
        {
            display_obj_mesh(tree3);
        }
    }
    glPopMatrix();
}


void draw_coin(float x, float y, float z)
{
   static GLfloat amt_dif[] = { 0.88, 0.85, 0.2};
   static GLfloat mat_specular[] = { 1.0, 0.1, 0.2, 1.0 };
   static GLfloat mat_shininess[] = { 150.0 };
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, amt_dif);
    
    glPushMatrix();
    {
        glTranslatef(x, y, z);
        
        
        glRotated(coin_angle, 0, 1, 0);
        glutSolidTorus(0.1, 0.5, 30, 30);
    }
    glPopMatrix();
}





void draw_ground_and_sky()
{
    glDisable(GL_LIGHTING);
    //draw the ground and street
    glEnable(GL_TEXTURE_2D);
    glPushMatrix();
    {
        int temp = user->z_pos;
        int p =  temp/view_boundry;
        
        for (int ii = 0; ii<p + 2; ii++)
        {
            const int subsection_num = 1;
            
            float bnd = view_boundry/subsection_num;// Draw Our Quad
            float sl = view_boundry;
            float ww = road_width;
            
            
            //sky
            glBindTexture(GL_TEXTURE_2D, sky_texture->texID);
            glBegin(GL_QUADS);
            glTexCoord2f(0.8f, 0.0f);    glVertex3f( ww*6,  street_height*3, -sl);
            glTexCoord2f(0.2f, 0.0f);    glVertex3f(-ww*6,  street_height*3, -sl);
            glTexCoord2f(0.2f, 1.0f);    glVertex3f(-ww*6,  street_height*3, 0);
            glTexCoord2f(0.8f, 1.0f);    glVertex3f( ww*6,  street_height*3, 0);
            glEnd();
            
            
            glPushMatrix();
            {
                
                for (int jj= 0; jj<subsection_num; jj++)
                {
                    
                    //ground
                    glBindTexture(GL_TEXTURE_2D, ground_texture->texID);
                    glBegin(GL_QUADS);
                    glTexCoord2f(0.8f, 0.0f);    glVertex3f( ww,  0.0f, -bnd);
                    glTexCoord2f(0.2f, 0.0f);    glVertex3f(-ww,  0.0f, -bnd);
                    glTexCoord2f(0.2f, 1.0f);    glVertex3f(-ww,  0.0f, 0);
                    glTexCoord2f(0.8f, 1.0f);    glVertex3f( ww,  0.0f, 0);
                    glEnd();
                    
                    glTranslatef(0, 0, -bnd);
                    
                }
            }glPopMatrix();
            
            glTranslatef(0, 0, -view_boundry);
        }
        
        // Done Drawing The Quad
    }glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    
    glEnable(GL_LIGHTING);
}

void draw_text()
{
    glLineWidth(4);
    glDisable(GL_LIGHTING);
    

    
    glPushMatrix();{
    if (start_counter>0)
    {
        char * str;
        if (start_counter>start_screen_time /4 *3)
        {
            str = "3";
        }
        else if (start_counter>start_screen_time/4*2 )
        {
            str = "2";
        }
        else if (start_counter > start_screen_time/4 )
        {
            str = "1";
        }else
        {
            str = "Go";
            millsecond_counter = 0;
        }
        
   
        glLineWidth(8);
        
        //DLOG("Let's get Ready");
        start_counter -= frame_time;
        glPushMatrix();
        {
            glColor3f(0.9, 0.3, 0.1);
            glTranslatef(-5.3,  8 +0.02, -user->z_pos - 10);
            
            glScalef(0.08, 0.1, 0.1);
            
            draw_str(str);
        }glPopMatrix();
    }
    else if (game_condition == GAME_ON)
    {
        glLineWidth(3);
        
        char dis_str[40];
        char str2[40];
        sprintf(dis_str, "Dist:%i/%i Spd:%i/%i Coin:%i", (int)user->z_pos, GOAL_DISTANCE,
                                                        (int)user->speed, (int)user->spd_max, user->coin);
        
        Player *fastestPc = engine->get_most_ahead_pc_player();
        
        if (user->z_pos >  fastestPc->z_pos)
        {
            sprintf(str2, "Time %i Ahead %i", get_elapsed_time(), (int)(user->z_pos - fastestPc->z_pos));
        }else
        {
            sprintf(str2, "Time %i Behind %i",  get_elapsed_time(), (int)(fastestPc->z_pos - user->z_pos));
        }
        
        const float font_x_offset = -6.5;
        float font_z_offset =  -user->z_pos - 9.2;
        
        
        glColor3f(0.8, 0.2, 0.3);
        //display user info
        glPushMatrix();
        {
            glTranslatef(font_x_offset,  11, font_z_offset);
            glScalef(0.004, 0.01, 0.01);
        
            draw_str(dis_str);
        }
        glPopMatrix();
        
        glPushMatrix();
        {
            glTranslatef(font_x_offset,  9, font_z_offset);
            glScalef(0.008, 0.01, 0.01);
            draw_str(str2);
            
        }glPopMatrix();
    }else
    {
        glLineWidth(10);
        
        char *str;
        if (game_condition == LOSE_THE_GAME)
        {
            str= "You lose....";
            
        }else
        {
            str = "You win!!";
        }
        
        
        
        glPushMatrix();
        {
            glTranslatef(-6,  9.2, -user->z_pos - 10);
            glScalef(0.02, 0.03, 0.031);
            glColor3f(0.8, 0.6, 0.28);
            draw_str(str);
            
        }glPopMatrix();
    }
    }glPopMatrix();
  
    glEnable(GL_LIGHTING);
}


void draw_instruction_page()
{
    
    glPushMatrix();
    {
        //printf("asds\n");
        //glcolor will affect texture
        glColor3f(1, 1, 1);
        glEnable(GL_TEXTURE_2D);
        
        glBindTexture(GL_TEXTURE_2D, instruction_page->texID);
        
        
        glBegin(GL_QUADS);
        glTexCoord2f(0.1f, 0.0f); glVertex2f(-road_width/2,  0);
        glTexCoord2f(0.1f, 0.9f); glVertex2f(-road_width/2,  street_height/2);
        glTexCoord2f(0.9f, 0.9f); glVertex2f( road_width/2,  street_height/2);
        glTexCoord2f(0.9f, 0.0f); glVertex2f( road_width/2,  0);
        glEnd();
        
        
        glDisable(GL_TEXTURE_2D);
        
    }glPopMatrix();
    
}

void draw_objects()
{
    
  
    
    draw_ground_and_sky();
    
    glEnable(GL_LIGHTING);
    
    
    //draw player
    glPushMatrix();
    {
        static bool draw = true;
        
        //flickering effect when hurt
        if (hurting_effect_counter >0)
        {
            draw = hurting_effect_counter%3 !=0;
            hurting_effect_counter--;
        }else
        {
            draw = true;
        }
        if (draw) {
            glTranslatef(get_opengl_x(user), 1+user->y_pos, get_opengl_z(user));
            
            glRotatef(-10, 1, 0, 0);
            glRotatef(90, 0, 1, 0);
            glScaled(0.5, 0.5, 0.5);
            glTranslatef(0, 2, 0);
            
            static float amt_color[] = {0.01, 0.01, 0.98};
            static float dif_color[] = {0.01, 0.01, 0.8};
            static float spc_color[] = {0.2, 0.0, 0.0};
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, dif_color);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, amt_color);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spc_color);
            
            humanList[0]->draw();
        }
    }
    glPopMatrix();
    
    
    for (int ii = 0; ii <engine->pc_player_list.size(); ii++)
    {
        Player *pc_player = engine->pc_player_list[ii];
        
        if (!isInUserViewRange(pc_player->z_pos))
        {
            continue;
        }
        
        //draw pc player
        glPushMatrix();
        {
            glTranslatef(get_opengl_x(pc_player), 1+pc_player->y_pos, get_opengl_z(pc_player));
            
            glRotatef(-10, 1, 0, 0);
            glRotatef(90, 0, 1, 0);
            glScaled(0.5, 0.5, 0.5);
            glTranslatef(0, 2, 0);
            
            float spc_color[] = {0.0, 0.4, 0.0};
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colors[ii]);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spc_color);
            
            humanList[ii+1]->draw();
            
        }
        glPopMatrix();
    }
    
   
    for (int ii = 0; ii<engine->pathObjectList.size(); ii++)
    {
        PathObject *obj = engine->pathObjectList.at(ii);
        
        //only draw visible object
        if (!isInUserViewRange(obj->z_pos))
        {
            continue;
        }
        
        //hit
        if (obj->ID == COIN_T)
        {
            draw_coin(get_opengl_x(obj), 1, get_opengl_z(obj));
        }
        else if ( obj->ID == SLOWDOWN_T)
        {
            draw_tree((SlowDownTrap *)obj);
        }
    }
    
    glDisable(GL_LIGHTING);
    //goal line
    if (isInUserViewRange(GOAL_DISTANCE))
    {
        const float scl = 0.03;
        
        glPushMatrix();{
        
        glColor3f(0.8, 0.1, 0.2);
        
        glTranslatef(0, 4.5, -GOAL_DISTANCE);
        glScalef(10, scl*2, scl);
        glutSolidCube(road_width+2);
        }glPopMatrix();
    }
    glEnable(GL_LIGHTING);
    
    
    //make the text in front of everything
    glClear(GL_DEPTH_BUFFER_BIT);
    draw_text();
    
    
    glColor3f(1, 1, 1);
    
}

void my_display()
{
    glClear(GL_COLOR_BUFFER_BIT |GL_DEPTH_BUFFER_BIT );
    
    coin_angle += 2;
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //DLOG(engine->player.x_pos);
    
    //make the camera a little bit behind the player
    my_cam->position.x = get_opengl_x(user);
    my_cam->position.z = get_opengl_z(user) + 3.3;
    my_cam->position.y = 3 + user->y_pos;
    
    my_cam->lookAt.x = my_cam->position.x;
    my_cam->lookAt.z = my_cam->position.z - view_boundry;
    my_cam->lookAt.y = 2;
    
    
    if (look_back)
    {
        my_cam->lookAt.z = -my_cam->lookAt.z;
    }
    
    my_cam->look();
    
    
    GLfloat light_amb_dif[] = {2.0, 1.0, 1.0, 1};
    
    GLfloat light_position1[] = { 0, 2, (user->z_pos) , 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
    glLightfv(GL_LIGHT0, GL_AMBIENT_AND_DIFFUSE, light_amb_dif);
    
    
    GLfloat light_position2[] = { -2, 1.0, (user->z_pos)+3, 0.0 };
    glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    
    
    
    
    draw_objects();
    glutSwapBuffers();
}


bool playing = false;
void *play_music_thread_method(void *filename)
{
    
    if (filename&&(!playing)&&music_on)
    {
        playing = true;
        
        char* name = (char *)filename;
        
#if defined(__APPLE__)
        char command[40];
        sprintf( command,"afplay %s", name);
        system(command);
#endif
        playing = false;
        
    }
    pthread_exit(NULL);
}


void play_music(char *fileName)
{
    pthread_t music_thread;
    pthread_create(&music_thread, NULL, play_music_thread_method, fileName);
    return;
}

void timeElapse(int value)
{
    int previous_coin = user->coin;
    float previous_spd = user->speed;
    
    
    
    if (!stop && start_counter <= 0) {
        
        millsecond_counter += frame_time;
        engine->move_forward(time_step);
        
        
        for (int ii = 0; ii < humanList.size(); ii++)
        {
            humanList[ii]->update_human_motion();
            humanList[ii]->update_human_motion();
            humanList[ii]->update_human_motion();
            humanList[ii]->update_human_motion();
            humanList[ii]->update_human_motion();
        }
        
        
        Player *best_pc_player = engine->get_most_ahead_pc_player();
        if (user->z_pos >= GOAL_DISTANCE&&best_pc_player->z_pos < GOAL_DISTANCE)
        {
            game_condition = WIN_THE_GAME;
            stop = true;
            printf("You win %i seconds\n", get_elapsed_time());
            
        }else if (user->z_pos < GOAL_DISTANCE&&best_pc_player->z_pos >= GOAL_DISTANCE)
        {
            game_condition = LOSE_THE_GAME;
            stop = true;
             printf("You lose %i seconds", get_elapsed_time());
        }else
        {
            game_condition = GAME_ON;
        }
        
    }
    

    if(user->coin > previous_coin)
    {
        play_music(coin_music_fn);
    }
    
    if (user->speed < previous_spd)
    {
        hurting_effect_counter = 60;
    }
    
       glutPostRedisplay();
    glutTimerFunc(frame_time, timeElapse, 0);
}


