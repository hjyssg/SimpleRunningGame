//
//  RunningGameEngine.cpp
//  Running
//
//  Created by Junyang Huang on 11/6/13.
//  Copyright (c) 2013 pitt. All rights reserved.
//

#include <vector>
#include  <set>
#include <iostream>

#include <string>
#include <fstream>
#include <sstream>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>

#if defined(_MSC_VER)
#include <windows.h>
#elif defined(__APPLE__)
#include <unistd.h>
#endif

#define  CLOG(a)  std::cout<<a<<"\n"
#if 0
#define  DLOG(a)  std::cout<<__PRETTY_FUNCTION__<<": "<<a<<"\n"
#else
#define DLOG(a)
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define COIN_T 1
#define COIN_T2 5
#define COIN_T3 6
#define SLOWDOWN_T 10
#define SLOWDOWN_T2 11
#define ATTACKER_T 100

#define DEFAULT_OBJECT_SIZE 1



const float top_speed_of_runner = 20;

//a parent class
class PathObject {
public:
    
    
    float z_pos;  //vertical coord
    int x_pos;  //horizontal coord
    float y_pos; //
    int radius;  //assume each object is a circle for collision detecting
    int ID;
    bool visible;
public:
    
 
    
    virtual std::string get_class_name()
    {
        return std::string("null object");
    }
    
    
    
};

class Chaser: public PathObject
{
public:
    float attack_point;
    float speed;
    
    Chaser()
    {
        radius = DEFAULT_OBJECT_SIZE;
        speed = 1;
        attack_point = 1;
    }
};

#define LEFT 1
#define RIGHT 0

//the moving object that user control
class Player: public PathObject
{
public:
    int HP;
    int coin;
    float speed; //only one parameter, becouse only move aloing -z axis
    float spd_max;
    float y_speed;
    float accel; //acceleration
    
    //these two toghter will decide how fast the ai to react the obstcale
    int AI_reaction_delay;
    int AI_reaction_ready;
    int AI_direction_preference;
    
    
public:
    Player()
    {
        radius = DEFAULT_OBJECT_SIZE;
        HP = 3;
        coin = 0;
        speed = 0;
        y_speed = 0;
        y_pos = 0;
        visible = true;
        spd_max = top_speed_of_runner;
        accel = 0.5;
        AI_reaction_ready = 0;
    }
    
    
    void move(float time_step)
    {
    	z_pos += speed * time_step;
        y_pos += y_speed * time_step;
        speed = MIN(speed+accel, spd_max);
        
        //jump
        if(y_pos>0)
        {
            y_speed = MIN(0, y_speed - 0.08);
        }else
        {
            y_speed = 0;
        }
    }
    
    //collect to increase spd
    //or increase max spd
    void spdUp_by_using_coin()
    {
        const int price = 3;
        
        if (coin < price)
        {
            return;
        }
        
    	if (speed < spd_max)
    	{
            speed = MIN(speed+1, spd_max);
            coin -= price;
    	}else
        {
            spd_max +=1;
            coin -= price;
        }
    }
    
    
};

class Coin: public PathObject
{
public:
    int value;
    
public:
    Coin(float h_p, float v_p)
    {
        radius = DEFAULT_OBJECT_SIZE;
        value = 1; //default is
        x_pos = h_p;
        z_pos = v_p;
        ID = COIN_T;
        visible = true;
    }
    
    std::string get_class_name()
    {
        return std::string("coin");
    }
    
};


#define SLOWDOWN_TYPE_1 0
#define SLOWDOWN_TYPE_2 1
#define SLOWDOWN_TYPE_3 2

class SlowDownTrap: public PathObject {
    
public:
    float friction;
    int type;
    
public:
    
    SlowDownTrap(float h_p, float v_p)
    {
        friction = 4; // default is 3
        radius = DEFAULT_OBJECT_SIZE;
        x_pos = h_p;
        z_pos = v_p;
        ID = SLOWDOWN_T;
        visible = true;
        type = rand()%3;
    }
    
    
    
    std::string get_class_name()
    {
        return std::string("slow-down trap");
    }
};




class Attacker: public PathObject
{
public:
    float attack_point;
    
public:
    Attacker(float h_p, float v_p)
    {
        attack_point = 1; //default
        radius = DEFAULT_OBJECT_SIZE;
        x_pos = h_p;
        z_pos = v_p;
        ID = ATTACKER_T;
        visible = true;
    }
    
    std::string get_class_name()
    {
        return std::string("attacker");
    }
};


class RunningWorld
{
    
public:
    int section_num;   //the num of section
    int section_range;  //the visible range
    int path_width; //the width of the path
    
    Player human_player;
    
    int difficuly;
    
    bool auto_mode; // auto_mode let ai to control the player runner
    
    
    std::vector<Player *>pc_player_list;
    
    
    std::vector<PathObject *> pathObjectList;
    
    RunningWorld(float sectionRange, int pathWidth, int ai_player_num)
    {
        section_range = sectionRange;
        section_num = 0;
        path_width = pathWidth;
        
        difficuly = 20;
        
        human_player = Player();
        human_player.x_pos = path_width/2.0;
        human_player.z_pos = 0;
        human_player.AI_reaction_delay = 4;
        human_player.AI_direction_preference = LEFT;
   
        auto_mode = false;
        
        for (int ii = 0; ii< ai_player_num; ii++)
        {
            Player *temp_pc_player = new Player();
            temp_pc_player->z_pos= 1;
            temp_pc_player->x_pos = path_width/4 * ii+1;
            temp_pc_player->ID = ii;
            temp_pc_player->AI_reaction_delay = 25;
            temp_pc_player->accel += 0.001*ii;
            
            if (rand()%2 == 0)
            {
                temp_pc_player->AI_direction_preference = RIGHT;
            }else
            {
                temp_pc_player->AI_direction_preference = LEFT;
            }
            
            
            pc_player_list.push_back(temp_pc_player);
            
        }
    }
    
    
    Player *get_most_behind_pc_player()
    {
        Player *slow_pc = pc_player_list[0];
        
        for (int ii = 1; ii < pc_player_list.size(); ii++)
        {
            if(pc_player_list[ii]->z_pos < slow_pc->z_pos)
            {
                slow_pc = pc_player_list[ii];
            }
        }
        
        return slow_pc;
    }
        
    Player *get_most_ahead_pc_player()
    {
        Player *fast_pc = pc_player_list[0];
        
        for (int ii = 1; ii < pc_player_list.size(); ii++)
        {
            if(pc_player_list[ii]->z_pos > fast_pc->z_pos)
            {
                fast_pc = pc_player_list[ii];
            }
        }
        
        return fast_pc;
    }
    
    void move_forward(float time_step)
    {
        human_player.move(time_step);
        for (int ii = 0; ii < pc_player_list.size(); ii++)
        {
            pc_player_list[ii]->move(time_step);
            
#if 1
            
            if (pc_player_list[ii]->AI_reaction_ready == pc_player_list[ii]->AI_reaction_delay)
            {
                decide_pc_player_next_move(pc_player_list[ii], false);
                pc_player_list[ii]->AI_reaction_ready = 0;
            }
            
            pc_player_list[ii]->AI_reaction_ready++;
            
#else
            decide_pc_player_next_move(pc_player_list[ii], false);
            
#endif
            
        }
        
        
        if (auto_mode)
        {
#if 1
     
            if (human_player.AI_reaction_ready == human_player.AI_reaction_delay)
            {
               decide_pc_player_next_move(&human_player, true);
                human_player.AI_reaction_ready = 0;
            }
            
            human_player.AI_reaction_ready++;
            
#else
            
            decide_pc_player_next_move(&human_player, true);

            
#endif
        }
        
        
        
        if (human_player.z_pos > section_num * section_range - section_range/3||get_most_ahead_pc_player()->z_pos > section_num * section_range- section_range/3)
        {
            for (int ii = 0; ii <section_range; ii++)
            {
                generate_random_objects(section_num * section_range + ii);
            }
            section_num ++;
        }
        
        erase_unvisible_object();
        detect_if_player_hit_anything(&human_player);
        for (int ii = 0; ii < pc_player_list.size(); ii++){
            detect_if_player_hit_anything(pc_player_list[ii]);
        }
        
    }
    
    //negative number move to the left
    //positive number move to the right
    bool move_player(Player *player_obj, int move)
    {
        int new_x = player_obj->x_pos +move;
        
        if (new_x >= 0 && new_x <= path_width)
        {
            player_obj->x_pos = new_x;
            return true;
        }
        
        return false;
    }
    
    
    //find the object neasest to the z_pos and along positive z axis
    PathObject * find_nearest_obj(float z_pos)
    {
        if (pathObjectList.size() ==0)
        {
            return NULL;
        }
        
        PathObject *nearest_obj = NULL;
        int ii = 0;
        for (ii = 0; ii < pathObjectList.size(); ii++)
        {
            //only stuff ahead of player
            if (pathObjectList[ii]->z_pos >= z_pos )
            {
                nearest_obj = pathObjectList[ii];
                break;
            }
        }
        
        if (nearest_obj != NULL)
        {
            for (int jj = ii; jj < pathObjectList.size(); jj++)
            {
                //only stuff ahead of player
                if (pathObjectList[jj]->z_pos >= z_pos && pathObjectList[jj]->z_pos < nearest_obj->z_pos)
                {
                    nearest_obj = pathObjectList[jj];
                }
            }
        }
        
        return nearest_obj;
    }
    
    //very AI for this game
    void decide_pc_player_next_move(Player *player_obj, bool for_user)
    {
        
        player_obj->spdUp_by_using_coin();
        
        PathObject *obj = find_nearest_obj(player_obj->z_pos);
        
        if (obj == NULL)
        {
            return;
        }
        
        if (obj->ID == COIN_T)
        {
            if (player_obj->x_pos < obj->x_pos )
            {
                //
                move_player(player_obj, 1);
                
            }else if (player_obj->x_pos > obj->x_pos )
            {
                move_player(player_obj, -1);
            }
            
        }else if (obj->ID == SLOWDOWN_T)
        {
            
            if (obj->x_pos == player_obj->x_pos)
            {
                
                if (for_user)
                {
                    DLOG("neasest is tree." << " user x, z"<<player_obj->x_pos << ","<<player_obj->z_pos <<".obj x, z" << obj->x_pos <<","<< obj->z_pos );
                }
                
                if (player_obj->AI_direction_preference == LEFT)
                {
                    
                    if (!move_player(player_obj, -1) )
                    {
                        move_player(player_obj, 1);
                    }
                }
                else
                {
                    if (!move_player(player_obj, 1) )
                    {
                        move_player(player_obj, -1);
                    }
                }
                
                
            }
        }
    }
    
    
    
    
    void player_jump(Player * player_obj)
    {
        player_obj->y_speed += 10;
    }
    
    
    void debug_print()
    {
        
        CLOG("human Player"<<" $"<<human_player.coin<<" (h,v) "<<human_player.x_pos<<","<<human_player.z_pos);
        
        for (int ii = 0; ii< pc_player_list.size(); ii++)
        {
            CLOG("Player"<<" $"<<pc_player_list[ii]->coin<<" (h,v) "<<pc_player_list[ii]->x_pos<<","<<pc_player_list[ii]->z_pos);
        }
        
        CLOG("====object list");
        for (int ii = 0; ii< pathObjectList.size(); ii++)
        {
            PathObject *obj = pathObjectList.at(ii);
            CLOG(obj->get_class_name()<<"(h,v) "<<obj->x_pos<<","<<obj->z_pos);
        }
    }
    
    
private:
    
    
    void detect_if_player_hit_anything(Player *player_obj)
    {
        if (pathObjectList.size()==0)
        {
            return;
        }
        
        for (std::vector<PathObject *>::iterator it = pathObjectList.begin() ; it != pathObjectList.end(); ++it)
        {
            PathObject *obj = *it;
            
            float h_distance = obj->z_pos - player_obj->z_pos;
             //no be hitted by object behind
            if (h_distance < 0)
            {
                continue;
            }
            
            
            if (obj->x_pos == player_obj->x_pos &&h_distance<(obj->radius+player_obj->radius))
            {
                //hit
                if (obj->ID == COIN_T)
                {
                    Coin * c = (Coin *)obj;
                    player_obj->coin += c->value;
                    
                    c->value = 0; //only collect gold once
                    c->visible = false;
                }
                else if ( obj->ID == SLOWDOWN_T)
                {
                    
                    player_obj->speed /= 10;
                    player_obj->z_pos = obj->z_pos+0.01;
                    
#if 0
                    SlowDownTrap *sl = (SlowDownTrap *) obj;
                    player_obj->speed -=sl->friction;
                    sl->friction = 0;
                    if (player_obj->speed<0) {
                        player_obj->speed = 0;
                    }  //can't move back
                    
#endif
                }
                else if (obj->ID == ATTACKER_T)
                {
                    Attacker *atk = (Attacker *) obj;
                    player_obj->HP -= atk->attack_point;
                    
                    atk->attack_point = 0; //only attack once
                }
                break;
            }
        }
    }
    
    void generate_random_objects(float z_position)
    {
        int token = rand()%difficuly;
        int h_c = rand()%(path_width+1);
        
        if (token == COIN_T||token == COIN_T2||token == COIN_T3)
        {
            PathObject * tempObj = new Coin(h_c, z_position);
            pathObjectList.push_back(tempObj);
        }
        else if (token ==SLOWDOWN_T||token ==SLOWDOWN_T2)
        {
            PathObject * tempObj = new SlowDownTrap(h_c, z_position);
            pathObjectList.push_back(tempObj);
        }
    }
    
    void erase_unvisible_object()
    {
        for (int ii = 0; ii< pathObjectList.size(); ii++)
        {
            PathObject *obj = pathObjectList.at(ii);
            
            if (obj->z_pos< MIN(human_player.z_pos-section_range, get_most_behind_pc_player()->z_pos - section_range) ||obj->visible == false)
            {
                pathObjectList.erase(pathObjectList.begin()+ii);
                delete obj;
            }
        }
    }
};





