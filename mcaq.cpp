/****************************************************************************

 //Adapted from example 5 of the Glui setup ---------------Pranav Dixit

****************************************************************************/


#include <algorithm>
#include <string.h>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <GL/glui.h>


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif




size_t f_id = 1;
size_t v_id = 1;

struct W_Edge;
struct W_Vertex;
struct W_Face;




struct W_Edge {


  W_Vertex *start, *end;
  W_Face *left, *right;
  W_Edge *left_prev, *left_next;
  W_Edge *right_prev, *right_next;
  W_Edge(){
  	start = nullptr;
  	left = nullptr;
    right = nullptr;
    left_prev = nullptr;
    left_next = nullptr;
    right_prev = nullptr;
    right_next = nullptr;
  }

};

 struct W_Vertex{

  int id;
  GLfloat x,y,z;
  W_Edge *edge;
  W_Vertex():id(1),x(0.0),y(0.0),z(0.0){

    edge = nullptr;
  }

};

struct W_Face{

  int id;
  W_Edge *edge;
  W_Face():id(1){
    edge = nullptr;
  }
};

//**********************************Class Declaration**************************************************

class Smf {

private:
  std::size_t num_vertices;
  std::size_t num_faces;

  std::map<size_t,W_Vertex*> vertices;
  std::map<size_t, W_Face* > faces;
  std::map<std::pair<size_t,size_t>, W_Edge* > edges;
  std::map<size_t, std::vector<GLfloat> > face_normals;
  std::map<size_t, std::vector<GLfloat> > vertex_normals;
  std::map<size_t, double > quadric_error_list;
  std::map<size_t, std::vector<std::vector<double> > > quadric_matrix_list;
  
  std::vector<W_Edge*> getEdges(const W_Vertex* vertex) const;
  std::vector<W_Face*> getFaces(const W_Vertex* vertex) const;
  std::vector<W_Edge*> getEdges(const W_Face* face) const ;
  std::vector<W_Vertex*> getVertices(const W_Face* face)const ;
 
//get face normals
  void getFaceNormal();
//get Vertex normals
  void getVertexNormal();

// get quadric list for all vertices
  bool getQuadricList(void);
// get quadric matrix for a edge
  std::vector<std::vector<double> > edgeQuadricMetric(W_Edge* edge);
// vertex quadric error for a vertex
  double getQuadricError( W_Vertex* vertex);
// get vertex quadric matrix for a vertex
  std::vector<std::vector<double> > getQuadricMatrix( W_Vertex* vertex);
// get new location for vertex
  std::vector<GLfloat> getNewLocation(const std::map<std::pair<size_t,size_t>, W_Edge* >::iterator &edge_it);
  // get matrix inverse
  std::vector<std::vector<double> > invertMatrix(const std::vector<std::vector<double> > &mat);


public:

// ~Smf(){}
friend std::ostream& operator<< (std::ostream& os, const Smf& smf);
~Smf(){}
Smf(const std::string &file= std::string());

void clearAll();
void conv_to_Winged(std::vector<size_t> face);
bool loadFile(const std::string &file);
bool saveFile(const std::string &file);
bool display();

int decimate(int k, int num_edges = 100);

};


//*************************************end of class declaration***********************************************

/** These are the  variables passed into GLUI ***************************************************/

float xy_aspect;
int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;

int   obj_type = 1;
int   light0_enabled = 1;
int   light1_enabled = 1;

int OPEN_FILE = 1;
int OUTPUT_FILE= 2;
int LOAD_MESH= 3;
int SAVE_FILE= 4;

char open_filetext[sizeof(GLUI_String)] = "cube";
char save_filetext[sizeof(GLUI_String)] = "test";

char open_filename[] = "./cube.smf";
char save_filename[] = "./test.smf";

char k[sizeof(GLUI_String)] = "10";
char collapse_edges[sizeof(GLUI_String)] = "1";

float light0_intensity = 1.0;
float light1_intensity = .4;
int   main_window;
float scale = 1.0;
int   show_mesh=1;
int   show_axes = 1;
int   show_text = 1;
float mesh_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

float obj_pos[] = { 0.0, 0.0, 0.0 };
const char *string_list[] = { "flat shaded", "smooth shaded", "wireframe", "shaded with mesh" };
int   curr_string = 3;

// ~Smf();

static Smf smf(open_filename);

/** Pointers to the windows and some of the controls we'll create **/
GLUI *glui, *glui2;
GLUI_Spinner    *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel;

/********** User IDs for callbacks ********/
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  260
#define ENABLE_ID            300
#define DISABLE_ID           301
#define SHOW_ID              302
#define HIDE_ID              303
#define SHADDING_ID          304
#define DECIMATE             305
#define K_VALUE              306
#define C_EDGES              307


/********** Miscellaneous global variables **********/

GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
GLfloat light0_position[] = {.5f, .5f, 1.0f, 0.0f};

GLfloat light1_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light1_diffuse[] =  {.9f, .6f, 0.0f, 1.0f};
GLfloat light1_position[] = {-1.0f, -1.0f, 1.0f, 0.0f};

GLfloat lights_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };




/*******************************************************************************************************************/


/*************************************Definition of class functions, constructor and friend functions*************/

//  Operator << is overloaded to check if model is uploaded properly

std::ostream& operator<< (std::ostream& os, const Smf& smf)
{


	for( std::map<size_t,W_Vertex* >::const_iterator it1 = smf.vertices.begin(); it1!= smf.vertices.end(); it1++){

      os << "v";

      os << " " << it1->second->x<<" "<< it1->second->y << " " << it1->second->z ;
     
      os <<std::endl; 


    }

    for( std::map<size_t, W_Face*>::const_iterator it2 = smf.faces.begin(); it2!= smf.faces.end(); it2++){

      os << "f"<< " " ;

      std::vector<W_Edge*> edges = smf.getEdges(it2->second);

      for (size_t i = 0; i < 3; i++) {
        if ((edges[i]->left) == (it2->second)) {
          os << edges[i]->start->id << " ";
        } else {
          os << edges[i]->end->id << " ";
        }
      }

      os<<std::endl;
    }

    // Uncomment below code, to see face normals and vertex normals

//   for(std::map<size_t,std::vector<GLfloat> >::const_iterator i = smf.vertex_normals.begin();i!= smf.vertex_normals.end(); i++){

//   os << "vn ";

//   for( std::vector<GLfloat>::const_iterator j  = i->second.begin(); j != i->second.end(); j++){

//     os << *j << " ";

//   }

//   os << std::endl;
// }

// for(std::map<size_t,std::vector<GLfloat>>::const_iterator i = smf.face_normals.begin();i!= smf.face_normals.end(); i++){

//   os << "fn";

//   for( std::vector<GLfloat>::const_iterator j  = i->second.begin(); j != i->second.end(); j++){

//     os << *j << " ";

//   }

//   os<<std::endl;
// }

for( std::map<std::pair<size_t,size_t>, W_Edge* >::const_iterator i = smf.edges.begin(); i != smf.edges.end(); i++){

  os << "e "<< i->second->start->id << " " << i->second->end->id<< std::endl;
}


return os;

}

/***start of function definitition from assignment1( loading SMF and coverting to Winged Edged data structure)**/



/************************************************************************************************************* 
SMF file handling and loading
*/
/***************************************************************************************************************/

//Constructor function

Smf::Smf(const std::string &file){

  if(file.length()==0){

    return;
  }
  // ~Smf();
  this->loadFile(file);

}


// clear all the data structures
void Smf::clearAll() {
  v_id = 1;
  f_id = 1;
  vertices.clear();
  faces.clear();
  edges.clear();
  face_normals.clear();
  vertex_normals.clear();
}

// member function to loadfile into model's data structures

bool Smf::loadFile(const std::string &file)

{

  std::ifstream ifile(file.c_str());

  std::cout<< "file: "<< file.c_str()<<std::endl;
  
  if(!ifile){
    std::cerr<< "Error occured while opening the file"<<std::endl;
    return false;
  }
  else{
    clearAll();
  }

std::string l;

while(std::getline(ifile, l)){

if(l.size()< 1 ){
  continue;
}

std::istringstream iss(l.substr(1));
std::vector<GLfloat> vertex(3,0.0);
std::vector<size_t> face(3,0);
std::vector<GLfloat> normal(3,0.0);

int check = 0;
switch(l[0]){

  case '#':

          iss>>num_vertices>>num_faces;

          // std::cout<< "Vertices: "<<num_vertices<<", Faces: "<<num_faces<<std::endl;

          break;

  case 'v':
          {

            iss>> vertex[0]>>vertex[1]>>vertex[2];
            W_Vertex* v = new W_Vertex;

            v->id = v_id;
            v->x = vertex[0];
            v->y = vertex[1];
            v->z = vertex[2];

            vertices.insert(std::make_pair(v_id,v));
            v_id++;
          
          }
          break;

  case 'f':
          {
          iss>> face[0]>>face[1] >>face[2];
          check++;
          conv_to_Winged(face);
          }
          break;
  default:

          break;

}


}

std::cout<< "Number of Vertices: "<<num_vertices<<" Check: "<<vertices.size()<<std::endl;
std::cout<<"Number of faces: "<<num_faces<<" Check: "<<faces.size()<<std::endl;
std::cout<<"Number of egdes: "<<(num_faces + num_vertices - 2)<<" Check : "<<edges.size()<<std::endl;

  if (vertices.size() != num_vertices || faces.size() != num_faces) 
  {
    std::cerr << "Check your file, if it is corrupted or not in the format specified in assignments" << file
        << std::endl;
    return false;
  }

  ifile.close();


getFaceNormal();
getVertexNormal();


std::cout<<"Face Normal Size: "<<face_normals.size()<<std::endl;
std::cout<<"Vertex Normal Size: "<<vertex_normals.size()<<std::endl;

std::cout<<"File: "<<file<< " loaded successfully"<<std::endl;

return true;

}

void Smf::conv_to_Winged(std::vector<size_t> face) {
 
  W_Face* w_face;
  w_face = new W_Face;
  W_Edge* w_edge[3]; //try to recreate the error you got when you used pointer to pointer to struct types, it was interesting to have that.

  w_edge[0] = new W_Edge;
  w_edge[1] = new W_Edge;
  w_edge[2] = new W_Edge;
  


  for (size_t i = 0; i < 2; ++i) {

    
    
    // checking if edge exist in edges.
    std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it = edges.find(std::make_pair(face[i + 1], face[i]));

    if (it == edges.end()) {

      w_edge[i]->start = (vertices.find(face[i])->second);

      w_edge[i]->end = (vertices.find(face[i + 1])->second);

      w_edge[i]->left = w_face;
      
      edges.insert(std::make_pair(std::make_pair(face[i], face[i + 1]), w_edge[i]));


      W_Vertex* w_vertex = vertices[face[(i + 1)]];

   
      if (w_vertex->edge == NULL){
      
        w_vertex->edge = w_edge[i];
      }

    } else {
      // if edge exists in edges then update right face of edge.

      it->second->right = w_face;
      w_edge[i] = it->second;
    }
  }

  // Handle back edge
  //checking if the edge exists
  std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it = edges.find(std::make_pair(face[0], face[2]));

  if (it == edges.end()) {
    w_edge[2]->start = (vertices.find(face[2])->second);
    w_edge[2]->end = (vertices.find(face[0])->second);
    w_edge[2]->left = w_face;
    edges.insert(std::make_pair(std::make_pair(face[2], face[0]), w_edge[2]));
    W_Vertex* w_vertex = vertices[face[0]];
    if (w_vertex->edge == NULL)
      w_vertex->edge = w_edge[2];
  } else {
    it->second->right = w_face;
    w_edge[2] = it->second;
  }


  // attaching edge to face
  w_face->edge = w_edge[0];
  w_face->id = f_id++;
  faces.insert(std::make_pair(w_face->id, w_face));

  // updating current edge's right_prev, right_next, left_prev and left_next
  for (size_t i = 1; i < 3; ++i) {
    std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it1 = edges.find(std::make_pair(face[i - 1], face[i]));
    std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it2 = edges.find(std::make_pair(face[i], face[i - 1]));
    if (it1 != edges.end()) {
      it1->second->left_prev = w_edge[i];
      it1->second->left_next = w_edge[(i + 1) % 3];
      //changing and checking
      // it1->second->left_prev = w_edge[(i + 1) % 3];
      // it1->second->left_next = w_edge[i];
    } else if (it2 != edges.end()) {
      it2->second->right_prev = w_edge[i];
      it2->second->right_next = w_edge[(i + 1) % 3];
      // it2->second->right_prev = w_edge[(i + 1) % 3];
      // it2->second->right_next = w_edge[i];
    }
  }
  // updating right_prev, right_next, left_prev and left_next for back edge
  std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it1 = edges.find(std::make_pair(face[2], face[0]));
  std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it2 = edges.find(std::make_pair(face[0], face[2]));
  if (it1 != edges.end()) {
    it1->second->left_prev = w_edge[0];
    it1->second->left_next = w_edge[1];
    // it1->second->left_prev = w_edge[1];
    // it1->second->left_next = w_edge[0];
  } else if (it2 != edges.end()) {
    it2->second->right_prev = w_edge[0];
    it2->second->right_next = w_edge[1];
    // it1->second->left_prev = w_edge[1];
    // it1->second->left_next = w_edge[0];
  
}  return ;
}




// saveFile: to save the model into a file
bool Smf::saveFile(const std::string &file){


  std::fstream f;

  f.open(file.c_str(), std::ios_base::out | std::ios_base::in);

  if(f.is_open())
  {
    std::cerr<< "File already exists, Please choose a different name";

    strcpy(save_filetext,"choose diff name");


    f.close();
  }
  else{

    f.clear();

    f.open(file.c_str(), std::ios_base::out);

    for( std::map<size_t,W_Vertex* >::iterator it1 = vertices.begin(); it1!= vertices.end(); it1++){

      f << "v";
    
      f << " " << it1->second->x<<" "<< it1->second->y << " " << it1->second->z ;

      f <<std::endl; 


    }

    for( std::map<size_t, W_Face*>::iterator it2 = faces.begin(); it2!= faces.end(); it2++){

      f << "f"<< " " ;

      std::vector<W_Edge*> edges = getEdges(it2->second);

      for (size_t i = 0; i < 3; i++) {
        if ((edges[i]->left) == (it2->second)) {
          f << edges[i]->start->id << " ";
        } else {
          f << edges[i]->end->id << " ";
        }
      }

      f<<std::endl;
    }

  }

  f.close();

  return true;
}

// get faces with given vertex
std::vector<W_Face*> Smf::getFaces(const W_Vertex* vertex) const {

  // std::cout<<"Vertex: "<<vertex<< ", vertex id:"<<vertex->id<<" : {} "<<vertex->x<<","<<vertex->y<<","<<vertex->y<<std::endl;
  W_Edge *e = vertex->edge;
  W_Edge *edge = e;
  // std::cout<<"vertex id: "<<vertex->id<<" ,vertex edge start id:"<< vertex->edge->start->id <<",vertex edge end id: "<<vertex->edge->end->id<<std::endl;
  std::vector<W_Face*> vFaces;
  do {

    if (edge->end == vertex) {
      vFaces.push_back(edge->right);
      edge = edge->right_next;
    } else {
      vFaces.push_back(edge->left);
      edge = edge->left_next;
    }
  } while (edge != e);
  return vFaces;
}


// get edges with given face
std::vector<W_Edge*> Smf::getEdges(const W_Face * face) const {

  W_Edge *e = face->edge;
  W_Edge *edge = e;
  std::vector<W_Edge*> vEdges;
  do {
    if (edge->right == face) {
      edge = edge->right_prev;
    } else {
      edge = edge->left_prev;
    }
    vEdges.push_back(edge);
  } while (edge != e);

  return vEdges;
}



// get vertices with given face
std::vector<W_Vertex*> Smf::getVertices(const W_Face* face) const {

   std::vector<W_Edge*> fEdges = getEdges(face);

    // std::cout<<"fEdges size: "<< fEdges.size()<<std::endl;
    std::vector<W_Vertex*> fVertices(3);
    for (int i = 0; i < 3; i++) 
    {
      if (fEdges[i]->left == face) {
        fVertices[i] = fEdges[i]->start;
      } else
        fVertices[i] = fEdges[i]->end;
    }

  return fVertices;
}



// get edges with given vertex
std::vector<W_Edge*> Smf::getEdges(const W_Vertex* vertex) const {
  W_Edge *e = vertex->edge;
  W_Edge *edge = e;
  std::vector<W_Edge*> vEdges;
  do {
    vEdges.push_back(edge);
    if (edge->end == vertex)
      edge = edge->right_next;
    else {
      edge = edge->left_next;
    }
  } while (!((edge->start == e->start && edge->end == e->end) || (edge->start == e->end && edge->end == e->start)));//(edge != e);
  return vEdges;
}



// get vertex normals
void Smf::getVertexNormal() {


  std::map<size_t, W_Vertex*>::iterator vertex_it = vertices.begin();
  
  size_t index = 0;
  while (vertex_it != vertices.end()) {
    index++;
    // std::cout<<"Vertex: "<<vertex_it->second<< " ,vertex id:"<<vertex_it->second->id<<" : {} "<<vertex_it->second->x<<","<<vertex_it->second->y<<","<<vertex_it->second->y<<std::endl;
    // std::cout<<"Vertex edge: "<< vertex_it->second->edge->start->id<<std::endl;

    std::vector<W_Face*> vFaces = getFaces(vertex_it->second);

    // std::cout<<"vertex index:"<<index<<"size"<<vFaces.size()<<std::endl;

    std::vector<GLfloat> vertexNormal(3, 0.0);

    for (size_t i = 0; i < vFaces.size(); i++) {
      for (size_t k = 0; k < 3; k++) {
        std::vector<GLfloat> temp = face_normals[vFaces[i]->id];
        vertexNormal[k] = vertexNormal[k] + /*temp[k];*/temp.at(k);
      }
    }

    // normalize it to one
    GLfloat length = sqrt(
        vertexNormal[0] * vertexNormal[0]
            + vertexNormal[1] * vertexNormal[1]
            + vertexNormal[2] * vertexNormal[2]);
    for (size_t k = 0; k < 3; k++) {
      vertexNormal[k] /= length;
    }
    vertex_normals.insert(make_pair(vertex_it->first, vertexNormal));
    vertex_it++;
  }
  return;
}


//get face normaals
void Smf::getFaceNormal() {
  std::map<size_t, W_Face*>::iterator face_it = faces.begin();
  while (face_it != faces.end()) {

  // std::cout<<"Are you here2"<<std::endl;
  std::vector<W_Edge*> fEdges = getEdges(face_it->second);

  // std::cout<<"fEdges size: "<< fEdges.size()<<std::endl;
  std::vector<W_Vertex*> fVertices(3);
  for (int i = 0; i < 3; i++) {
    if (fEdges[i]->left == face_it->second) {
      fVertices[i] = fEdges[i]->start;
    } else
      fVertices[i] = fEdges[i]->end;
  }


  std::vector<GLfloat> diff01(3,0);
  std::vector<GLfloat> diff12(3,0);
  std::vector<GLfloat> normal(3,0);
  GLfloat length;




  diff01[0] = fVertices[1]->x - fVertices[0]->x;
  diff12[0] = fVertices[2]->x - fVertices[1]->x;

  diff01[1] = fVertices[1]->y - fVertices[0]->y;
  diff12[1] = fVertices[2]->y - fVertices[1]->y;

  diff01[2] = fVertices[1]->z - fVertices[0]->z;
  diff12[2] = fVertices[2]->z - fVertices[1]->z;


  normal[0] = diff01[1]*diff12[2]- diff01[2]*diff12[1];
  normal[1] = diff01[2]*diff12[0]- diff01[0]*diff12[2];
  normal[2] = diff01[0]*diff12[1]- diff01[1]*diff12[0];


   length = std::sqrt(std::pow(normal[0],2)+std::pow(normal[1],2)+std::pow(normal[2],2));

  for (int i = 0 ; i < 3; i++){
    normal[i]/=length;
  }

  face_normals.insert(make_pair(face_it->first, normal));
  face_it++;
  }
  
  return;
}




bool Smf::display(){

glBegin(GL_TRIANGLES);

std::vector<GLfloat> normal(3);
std::vector<W_Edge*> fEdges(3);

 for(std::map<size_t, W_Face*>::iterator face_it = faces.begin();face_it!=faces.end();++face_it){ 

      W_Face* f = face_it->second;

      fEdges = getEdges(f);

      for( size_t j = 0 ; j < fEdges.size(); j++){


        switch(curr_string){
    // handling flat faced model, flat models using face normal even when vertex is shared across faces.
          case 0: 
          normal = face_normals[face_it->first];
          break;

          default: 

          if(fEdges[j]->left==f)
            normal = vertex_normals[fEdges[j]->start->id];
          else 
            normal = vertex_normals[fEdges[j]->end->id];

          break;
        }

      //Normalizing 

      GLfloat length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
      for (size_t k = 0; k < 3; ++ k)
      {
        normal[k] /= length;
      }
        
      glNormal3f(normal[0],normal[1],normal[2]);

      

      if(fEdges[j]->left==f){

        glVertex3f(fEdges[j]->start->x, fEdges[j]->start->y,fEdges[j]->start->z);
        
        } else {

          glVertex3f(fEdges[j]->end->x, fEdges[j]->end->y, fEdges[j]->end->z);      
      }

        

      }

      
    }

glEnd();

return true;

}


/********************************************************************************************************
/* End of Smf class definitions and functions


**********************************************************************************************************/


/*************Decimattion functions:******************/

/******************************************************************************************************/



int Smf::decimate(int k, int num_edges)
{
if( edges.size() < num_edges){

  std::cout<<"Number of edges are less than the specified number to be removed"<<std::endl;
  return EXIT_FAILURE;
}

  // get quadric error and matrix list for each vertex
  getQuadricList();

  std::cout<<"successfully created quadric list"<<std::endl;

  std::cout<<"Matrix list size: "<<quadric_matrix_list.size()<<std::endl;
  std::cout<<"Error list size: "<< quadric_error_list.size()<<std::endl;
  
  size_t original_size = edges.size();
  srand( time(NULL) );

  // run while loop till specified number of edges(num_edges) are removed

  while (edges.size() > (original_size - num_edges))
  {
    std::set<size_t> candidates;

    // choose k candidate edges randomly
    for (size_t i = 0; i < k; ++ i)
    {
      size_t r = std::rand() % edges.size();
      candidates.insert(r);
    }


    double least_quadric_error = DBL_MAX;
    double prev_least_error = DBL_MAX;
    std::set<size_t>::iterator prev_least_it;
    std::set<size_t>::iterator least;

    // find the candidate edge with least error
    for (std::set<size_t>::iterator it = candidates.begin(); it != candidates.end(); ++ it)
    {
      // int flag = 0;
      std::map<std::pair<size_t,size_t>, W_Edge* >::iterator edge_it = edges.begin();
      std::advance(edge_it, *it);

      std::vector<W_Edge*> start_Edges = getEdges(edge_it->second->start);
      std::vector<W_Edge*> end_Edges = getEdges(edge_it->second->end);
      if(start_Edges.size()<4 || end_Edges.size()<4 ||start_Edges.size()>5 || end_Edges.size()>5){

      // std::cout<<"ignoring the decimation of a edge with points having valence less than 4"<<std::endl;
      continue;

    }
      double quadric = quadric_error_list[ edge_it->second->start->id ] + quadric_error_list[ edge_it->second->end->id];
      if (quadric < least_quadric_error)
      {
        // prev_least_error = least_quadric_error;
        least_quadric_error = quadric;
        // prev_least_it = least;
        least = it;
        // flag =1;
      }


    }

   


// point to the edge with least error (least_edge->second)
    std::map<std::pair<size_t,size_t>, W_Edge* >::iterator least_edge = edges.begin();
    std::advance(least_edge, *least);

      std::vector<W_Edge*> start_Edges = getEdges(least_edge->second->start);
      std::vector<W_Edge*> end_Edges = getEdges(least_edge->second->end);

      if(start_Edges.size()<4 || end_Edges.size()<4 ||start_Edges.size()>5 || end_Edges.size()>5){

      // std::cout<<"ignoring the decimation of a edge with points having valence less than 4"<<std::endl;
      continue;

    }

// get the new location for the start vertex
   
    std::vector<GLfloat> newVertex = getNewLocation(least_edge);
    least_edge->second->start->x = newVertex[0];
    least_edge->second->start->y = newVertex[1];
    least_edge->second->start->z = newVertex[2];
    if (least_edge->second->start->edge == least_edge->second) {
      least_edge->second->start->edge = least_edge->second->left_next;
    }

    // std::cout<<smf<<std::endl;

    std::cout<<"Deleting Edge: "<< least_edge->second->start->id<<"-----"<<least_edge->second->end->id<<std::endl;

//update faces and Vertices involved in the decimation

  std::vector<W_Edge*> vEdges = getEdges(least_edge->second->end);
  for (int i = 0; i < vEdges.size(); i++) {

    std::map<std::pair<size_t, size_t>, W_Edge*>::iterator it1 = edges.find(
                                                      std::make_pair(vEdges[i]->start->id, vEdges[i]->end->id));


// ignoring the least edge as it will be deleted
    if (vEdges[i]->start == least_edge->second->start && vEdges[i]->end == least_edge->second->end) {
      continue;
    } 
    else if (least_edge->second->left_prev->start == vEdges[i]->start 
                                                      && least_edge->second->left_prev->end == vEdges[i]->end) {

     
      if (least_edge->second->left_prev->left->edge == least_edge->second->left_prev) {
         // update face edge
        least_edge->second->left_prev->left->edge = least_edge->second->left_prev->left_prev;

      // update vertex edge
        if (least_edge->second->left_prev->end->edge == least_edge->second->left_prev) {
          least_edge->second->left_prev->end->edge = least_edge->second->left_prev->left_prev;
        }

      }
      
      if (least_edge->second->left_prev->right->edge == least_edge->second->left_prev) {
        // update face edge
        least_edge->second->left_prev->right->edge = least_edge->second->left_prev->right_prev;

        // update vertex edge
        if (least_edge->second->left_prev->start->edge == least_edge->second->left_prev) {
          least_edge->second->left_prev->start->edge = least_edge->second->left_prev->right_prev;
        }
      }

      continue;

    } 
    else if (least_edge->second->right_next->start == vEdges[i]->start
                                                        && least_edge->second->right_next->end == vEdges[i]->end) {
      
      if (least_edge->second->right_next->left->edge == least_edge->second->right_next) {
        // update face edge
        least_edge->second->right_next->left->edge = least_edge->second->right_next->left_next;

        // update vertex edge
        if (least_edge->second->right_next->start->edge == least_edge->second->right_next) {
          least_edge->second->right_next->start->edge = least_edge->second->left_next;
        }

      }

      if (least_edge->second->right_next->right->edge == least_edge->second->right_next) {
        // update face edge
        least_edge->second->right_next->right->edge = least_edge->second->right_next->right_next;
        // update vertex edge
        if (least_edge->second->right_next->end->edge == least_edge->second->right_next) {
          least_edge->second->right_next->end->edge = least_edge->second->right_next;
        }
      }

      continue;
    } 
    else if (vEdges[i]->start == least_edge->second->end) {

  // directing the start of the edge to least_edge->second-start
      it1->second->start = least_edge->second->start;

    } else if (vEdges[i]->end == least_edge->second->end) {
  // directing the end of the edge to least_edge->second-start
      it1->second->end = least_edge->second->start;
    }

     // remove edge
    edges.erase(it1->first);//check map erase

    // inserting updated edge as the pair needs to change
    edges.insert(
        std::make_pair(
            std::make_pair(it1->second->start->id, it1->second->end->id),
            it1->second));


  }


    // remove faces
    faces.erase(least_edge->second->left->id);
    faces.erase(least_edge->second->right->id);

    // remove edges
    edges.erase(std::make_pair(least_edge->second->start->id, least_edge->second->end->id));
    edges.erase(std::make_pair(least_edge->second->left_prev->start->id, least_edge->second->left_prev->end->id));
    edges.erase(std::make_pair(least_edge->second->right_next->start->id,least_edge->second->right_next->end->id));

    // remove vertex
    vertices.erase(least_edge->second->end->id);


    //Updating Edges************************************************************

    //left hand side edges
    W_Edge* l0 = least_edge->second->left_prev;
    W_Edge* l1 = least_edge->second->left_next;

//setting directions of the edges
  if (l1->right == least_edge->second->left) {
    if (l0->left == least_edge->second->left) {
      l1->right = l0->right;
      l1->right_next = l0->right_next;
      l1->right_prev = l0->right_prev;
  
      if (l0->right_next->right == l0->right) {
        l0->right_next->right_prev = l1;
      } else {
        l0->right_next->left_prev = l1;
      }
      if (l0->right_prev->right == l0->right) {
        l0->right_prev->right_next = l1;
      } else {
        l0->right_prev->left_next = l1;
      }

    } else {
      // if l0's left face is not same as least edge's left
      l1->right = l0->left;
      l1->right_next = l0->left_next;
      l1->right_prev = l0->left_prev;

      if (l0->left_next->left == l0->left) {
        l0->left_next->left_prev = l1;
      } else {
        l0->left_next->right_prev = l1;
      }
      if (l0->left_prev->left == l0->left) {
        l0->left_prev->left_next = l1;
      } else {
        l0->left_prev->right_next = l1;
      }
    }
  } // taking opposite direction of l1
  else if (l1->left == least_edge->second->left) {
    if (l0->left == least_edge->second->left) {
      l1->left = l0->right;
      l1->left_next = l0->right_next;
      l1->left_prev = l0->right_prev;

      if (l0->right_next->right == l0->right) {
        l0->right_next->right_prev = l1;
      } else {
        l0->right_next->left_prev = l1;
      }
      if (l0->right_prev->right == l0->right) {
        l0->right_prev->right_next = l1;
      } else {
        l0->right_prev->left_next = l1;
      }

    } else {
      l1->left = l0->left;
      l1->left_next = l0->left_next;
      l1->left_prev = l0->left_prev;


      if (l0->left_next->left == l0->left) {
        l0->left_next->left_prev = l1;
      } else {
        l0->left_next->right_prev = l1;
      }
      if (l0->left_prev->left == l0->left) {
        l0->left_prev->left_next = l1;
      } else {
        l0->left_prev->right_next = l1;
      }

    }

  }

  // right hand side edges
  W_Edge* r0 = least_edge->second->right_prev;
  W_Edge* r1 = least_edge->second->right_next;

// setting directions of the edges
  // if r0's left is least edge's right
  if (r0->left == least_edge->second->right) {
    if (r1->left == least_edge->second->right) {
      r0->left = r1->right;
      r0->left_prev = r1->right_prev;
      r0->left_next = r1->right_next;

      if (r1->right_next->right == r1->right) {
        r1->right_next->right_prev = r0;
      } else {
        r1->right_next->left_prev = r0;
      }
      if (r1->right_prev->right == r1->right) {
        r1->right_prev->right_next = r0;
      } else {
        r1->right_prev->left_next = r0;
      }
    } else {
      r0->left = r1->left;
      r0->left_prev = r1->left_prev;
      r0->left_next = r1->left_next;

    
      if (r1->left_next->left == r1->left) {
        r1->left_next->left_prev = r0;
      } else {
        r1->left_next->right_prev = r0;
      }
      if (r1->left_prev->left == l0->left) {
        r1->left_prev->left_next = r0;
      } else {
        r1->left_prev->right_next = r0;
      }
    }

//if r0's right is least edge's right face
  } else if (r0->right == least_edge->second->right) {
    // if r1's left is least edge's right face
    if (r1->left == least_edge->second->right) {
      r0->right = r1->right;
      r0->right_prev = r1->right_prev;
      r0->right_next = r1->right_next;

      // updating nei
      if (r1->right_next->right == r1->right) {
        r1->right_next->right_prev = r0;
      } else {
        r1->right_next->left_prev = r0;
      }
      if (r1->right_prev->right == r1->right) {
        r1->right_prev->right_next = r0;
      } else {
        r1->right_prev->left_next = r0;
      }

    }
    // if r1's right is least edge's right face 
    else {
      r0->right = r1->left;
      r0->right_prev = r1->left_prev;
      r0->right_next = r1->left_next;
      // updating neighbours
      if (r1->left_next->left == r1->left) {
        r1->left_next->left_prev = r0;
      } else {
        r1->left_next->right_prev = r0;
      }
      if (r1->left_prev->left == l0->left) {
        r1->left_prev->left_next = r0;
      } else {
        r1->left_prev->right_next = r0;
      }

    }
  }

  
//update  quadrics

    quadric_matrix_list[least_edge->second->start->id] = edgeQuadricMetric(least_edge->second);
    quadric_error_list[least_edge->second->start->id] = least_quadric_error;
    quadric_matrix_list.erase(least_edge->second->end->id);
    quadric_error_list.erase(least_edge->second->end->id);


// std::cout<<smf<<std::endl;

}

  
  return EXIT_SUCCESS;
}




// get quadric error and matrix list for each vertex
bool Smf::getQuadricList(void)
{
  quadric_matrix_list.clear();
  quadric_error_list.clear();
  // for each vertex
  
  for (std::map<size_t, W_Vertex*>::iterator it = vertices.begin();it!=vertices.end();++it)
  {
    

    std::vector<std::vector<double> > quad = getQuadricMatrix(it->second);

    // std::cout<<"quad:[33]: "<< quad[3][3]<<std::endl;
    quadric_matrix_list.insert(std::make_pair(it->first,quad ));

    // std::cout<< "matrix list size: "<<quadric_matrix_list.size()<<std::endl;

    quadric_error_list.insert( std::make_pair(it->first, getQuadricError(it->second)) );
  }
  return true;
}

std::vector<std::vector<double> > Smf::edgeQuadricMetric(W_Edge* edge) {
 

  std::vector<std::vector<double> > U = quadric_matrix_list[edge->start->id];
  std::vector<std::vector<double> > W = quadric_matrix_list[edge->end->id];


  std::vector<std::vector<double> > Q(4, std::vector<double>(4, 0.0));

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Q[i][j] = U[i][j] + W[i][j];
    }
  }
  Q[3][0] = Q[3][1] = Q[3][2] = 0.0;
  Q[3][3] = 1.0;
  
  return Q;

}



// vertex quadric error
double Smf::getQuadricError( W_Vertex* vertex)
{
 

  std::map<size_t, W_Vertex*>::iterator it = vertices.find(vertex->id);
    if (it == vertices.end()){
    return DBL_MAX;

    }

  // get quadric error matrix
  std::vector<std::vector<double> > Q = quadric_matrix_list[vertex->id];

  // get x,y,z of the vertex in 4 length vector
  std::vector<GLfloat> v = {vertex->x,vertex->y,vertex->z,1.0};
  

  double error = 0.0;
  // error
  for (size_t i = 0; i < 4; ++ i)
  {
    double e = 0.0;
    for (size_t j = 0; j < 4; ++ j)
    {
      e += v[j] * Q[j][i];
    }
    error += e * v[i];
  }
  // std::cout<<"Error:"<<error<<std::endl;
  return error;
}


// get vertex quadric matrix
std::vector<std::vector<double> > Smf::getQuadricMatrix( W_Vertex* vertex)

{

  std::map<size_t, W_Vertex*>::iterator it = vertices.find(vertex->id);
    if (it == vertices.end()){
    return std::vector<std::vector<double> >();

    }

  // quadric matrix
  std::vector<std::vector<double> > Q(4, std::vector<double> (4, 0.0));
  double plane[4] = {0.0, 0.0, 0.0, 0.0};


  // get faces for a vertex

  // get all the vertexes of the face

  // get the equation of the plane of that face with vertexes

  // calculate Q for the equation of the plane or from a, b, c, d

  // for connected faces

  std::vector<W_Face*> vertex_faces = getFaces(vertex);


  for(int i = 0; i<vertex_faces.size();i++)
  {

    std::vector<std::vector<double> > Kp(4, std::vector<double> (4, 0.0));

    std::vector<W_Vertex*>  face_vertices = getVertices(vertex_faces[i]);

    // std::cout<<"Face : "<< vertex_faces[i]->id<<" , "<< "Vertices : "<<face_vertices[0]->id<<", "<<face_vertices[1]->id<<", "<<face_vertices[2]->id<<std::endl;

   
    if(face_vertices.size()!=3){

      std::cout<<"Something is wrong with getVertices function"<<std::endl;
    }

    std::vector<GLfloat> A = {face_vertices[0]->x, face_vertices[1]->x, face_vertices[2]->x};

    std::vector<GLfloat> B = {face_vertices[0]->y, face_vertices[1]->y, face_vertices[2]->y};

    std::vector<GLfloat> C = {face_vertices[0]->z, face_vertices[1]->z, face_vertices[2]->z};

    // calculate parameters of a plane (a, b, c, d)
    plane[0] = (B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]);
    plane[1] = (B[2] - A[2]) * (C[0] - A[0]) - (C[2] - A[2]) * (B[0] - A[0]);
    plane[2] = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);


    // normalize
    double normalizer = sqrt( plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2] );
    plane[0] /= normalizer;
    plane[1] /= normalizer;
    plane[2] /= normalizer;
    plane[3] = - (plane[0] * A[0] + plane[1] * A[1] + plane[2] * A[2]);
    for (size_t i = 0; i < 4; ++ i)
    {
      for (size_t j = 0; j < 4; ++ j)
      {
   // calculate Kp
        Kp[i][j] = plane[i] * plane[j];
        // sum to get quadric
        Q[i][j] += Kp[i][j];

      
      }
    }
  }

  return Q;
}


// get new location for vertex
std::vector<GLfloat> Smf::getNewLocation(const std::map<std::pair<size_t,size_t>, W_Edge* >::iterator &edge_it)
{
  

  std::vector<std::vector<double> > Q = edgeQuadricMetric(edge_it->second);
  // invert the matrix
  std::vector<std::vector<double> > inv =  invertMatrix(Q);

  std::vector<GLfloat> new_vertex(3, 0.0);

  // not invertible
  if (inv.size() < 4)
  {
  
    new_vertex[0] = (vertices[edge_it->second->start->id]->x + vertices[edge_it->second->end->id]->x) / 2.0;
    new_vertex[1] = (vertices[edge_it->second->start->id]->y + vertices[edge_it->second->end->id]->y) / 2.0;
    new_vertex[2] = (vertices[edge_it->second->start->id]->z + vertices[edge_it->second->end->id]->z) / 2.0;
    
    // new_vertex[3] = 1.0;
  }
  else
  {
    // if inverse exist
    
    new_vertex[0] = inv[0][3];
    new_vertex[1] = inv[1][3];
    new_vertex[2] = inv[2][3];
  }

  return new_vertex;
}


//found the below code from external resource to calculate the inverse
std::vector<std::vector<double> > Smf::invertMatrix(const std::vector<std::vector<double> > &mat) {
 

  std::vector<std::vector<double> > inverse(mat.size(),
      std::vector<double>(mat.size(), 0.0));
  double detrminant =0.0;
// for 4 x 4 mat
  inverse[0][0] = mat[1][1] * mat[2][2] * mat[3][3]
      + mat[1][2] * mat[2][3] * mat[3][1]
      + mat[1][3] * mat[2][1] * mat[3][2]
      - mat[1][1] * mat[2][3] * mat[3][2]
      - mat[1][2] * mat[2][1] * mat[3][3]
      - mat[1][3] * mat[2][2] * mat[3][1];
  inverse[0][1] = mat[0][1] * mat[2][3] * mat[3][2]
      + mat[0][2] * mat[2][1] * mat[3][3]
      + mat[0][3] * mat[2][2] * mat[3][1]
      - mat[0][1] * mat[2][2] * mat[3][3]
      - mat[0][2] * mat[2][3] * mat[3][1]
      - mat[0][3] * mat[2][1] * mat[3][2];
  inverse[0][2] = mat[0][1] * mat[1][2] * mat[3][3]
      + mat[0][2] * mat[1][3] * mat[3][2]
      + mat[0][3] * mat[1][1] * mat[3][2]
      - mat[0][1] * mat[1][3] * mat[3][2]
      - mat[0][2] * mat[1][1] * mat[3][3]
      - mat[0][3] * mat[1][2] * mat[3][1];
  inverse[0][3] = mat[0][1] * mat[1][3] * mat[2][2]
      + mat[0][2] * mat[1][1] * mat[2][3]
      + mat[0][3] * mat[1][2] * mat[2][1]
      - mat[0][1] * mat[1][2] * mat[2][3]
      - mat[0][2] * mat[1][3] * mat[2][1]
      - mat[0][3] * mat[1][1] * mat[2][2];
  inverse[1][0] = mat[1][0] * mat[2][3] * mat[3][2]
      + mat[1][2] * mat[2][1] * mat[3][3]
      + mat[1][3] * mat[2][2] * mat[3][0]
      - mat[1][0] * mat[2][2] * mat[3][3]
      - mat[1][2] * mat[2][3] * mat[3][0]
      - mat[1][3] * mat[2][0] * mat[3][2];
  inverse[1][1] = mat[0][0] * mat[2][2] * mat[3][3]
      + mat[0][2] * mat[2][3] * mat[3][0]
      + mat[0][3] * mat[2][0] * mat[3][2]
      - mat[0][0] * mat[2][3] * mat[3][2]
      - mat[0][2] * mat[2][0] * mat[3][3]
      - mat[0][3] * mat[2][2] * mat[3][0];
  inverse[1][2] = mat[0][0] * mat[1][3] * mat[3][2]
      + mat[0][2] * mat[1][0] * mat[3][3]
      + mat[0][3] * mat[1][2] * mat[3][0]
      - mat[0][0] * mat[1][2] * mat[3][3]
      - mat[0][2] * mat[1][3] * mat[3][0]
      - mat[0][3] * mat[1][0] * mat[3][2];
  inverse[1][3] = mat[0][0] * mat[1][2] * mat[2][3]
      + mat[0][2] * mat[1][3] * mat[2][0]
      + mat[0][3] * mat[1][0] * mat[2][2]
      - mat[0][0] * mat[1][3] * mat[2][2]
      - mat[0][2] * mat[1][0] * mat[2][3]
      - mat[0][3] * mat[1][2] * mat[2][0];
  inverse[2][0] = mat[1][0] * mat[2][1] * mat[3][3]
      + mat[1][1] * mat[2][3] * mat[3][0]
      + mat[1][3] * mat[2][0] * mat[3][1]
      - mat[1][0] * mat[2][3] * mat[3][1]
      - mat[1][1] * mat[2][0] * mat[3][3]
      - mat[1][3] * mat[2][1] * mat[3][0];
  inverse[2][1] = mat[0][0] * mat[2][3] * mat[3][1]
      + mat[0][1] * mat[2][0] * mat[3][3]
      + mat[0][3] * mat[2][1] * mat[3][0]
      - mat[0][0] * mat[2][1] * mat[3][3]
      - mat[0][1] * mat[2][3] * mat[3][0]
      - mat[0][3] * mat[2][0] * mat[3][1];
  inverse[2][2] = mat[0][0] * mat[1][1] * mat[3][3]
      + mat[0][1] * mat[1][3] * mat[3][0]
      + mat[0][3] * mat[1][0] * mat[3][1]
      - mat[0][0] * mat[0][3] * mat[3][1]
      - mat[0][1] * mat[1][0] * mat[3][3]
      - mat[0][3] * mat[1][1] * mat[3][0];
  inverse[2][3] = mat[0][0] * mat[1][3] * mat[2][1]
      + mat[0][1] * mat[1][0] * mat[2][3]
      + mat[0][3] * mat[1][1] * mat[2][0]
      - mat[0][0] * mat[1][1] * mat[2][3]
      - mat[0][1] * mat[1][3] * mat[2][0]
      - mat[0][3] * mat[1][0] * mat[2][1];
  inverse[3][0] = mat[1][0] * mat[2][2] * mat[3][1]
      + mat[1][1] * mat[2][0] * mat[3][2]
      + mat[1][2] * mat[2][1] * mat[3][0]
      - mat[1][0] * mat[2][1] * mat[3][2]
      - mat[1][1] * mat[2][2] * mat[3][0]
      - mat[1][2] * mat[2][0] * mat[3][1];
  inverse[3][1] = mat[0][0] * mat[2][1] * mat[3][2]
      + mat[0][1] * mat[2][2] * mat[3][0]
      + mat[0][2] * mat[2][0] * mat[3][1]
      - mat[0][0] * mat[2][2] * mat[3][1]
      - mat[0][1] * mat[2][0] * mat[3][2]
      - mat[0][2] * mat[2][1] * mat[3][0];
  inverse[3][2] = mat[0][0] * mat[1][2] * mat[3][1]
      + mat[0][1] * mat[1][0] * mat[3][2]
      + mat[0][2] * mat[1][1] * mat[3][0]
      - mat[0][0] * mat[1][1] * mat[3][2]
      - mat[0][1] * mat[1][2] * mat[3][0]
      - mat[0][2] * mat[1][0] * mat[3][1];
  inverse[3][3] = mat[0][0] * mat[1][1] * mat[2][2]
      + mat[0][1] * mat[1][2] * mat[2][0]
      + mat[0][2] * mat[1][0] * mat[2][1]
      - mat[0][0] * mat[1][2] * mat[2][1]
      - mat[0][1] * mat[1][0] * mat[2][2]
      - mat[0][2] * mat[1][1] * mat[2][0];
// get determinant
  double determinant = mat[0][0] * mat[1][1] * mat[2][2]
      * mat[3][3]
      + mat[0][0] * mat[1][2] * mat[2][3] * mat[3][1]
      + mat[0][0] * mat[1][3] * mat[2][1] * mat[3][2] +

      mat[0][1] * mat[1][0] * mat[2][3] * mat[3][2]
      + mat[0][1] * mat[1][2] * mat[2][0] * mat[3][3]
      + mat[0][1] * mat[1][3] * mat[2][2] * mat[3][0] +

      mat[0][2] * mat[1][0] * mat[2][1] * mat[3][3]
      + mat[0][2] * mat[1][1] * mat[2][3] * mat[3][0]
      + mat[0][2] * mat[1][3] * mat[2][0] * mat[3][1] +

      mat[0][3] * mat[1][0] * mat[2][2] * mat[3][1]
      + mat[0][3] * mat[1][1] * mat[2][0] * mat[3][2]
      + mat[0][3] * mat[1][2] * mat[2][1] * mat[3][0] -

  mat[0][0] * mat[1][1] * mat[2][3] * mat[3][2]
      - mat[0][0] * mat[1][2] * mat[2][1] * mat[3][3]
      - mat[0][0] * mat[1][3] * mat[2][2] * mat[3][1] -

      mat[0][1] * mat[1][0] * mat[2][2] * mat[3][3]
      - mat[0][1] * mat[1][2] * mat[2][3] * mat[3][0]
      - mat[0][1] * mat[1][3] * mat[2][0] * mat[3][2] -

      mat[0][2] * mat[1][0] * mat[2][3] * mat[3][1]
      - mat[0][2] * mat[1][1] * mat[2][0] * mat[3][3]
      - mat[0][2] * mat[1][3] * mat[2][1] * mat[3][0] -

      mat[0][3] * mat[1][0] * mat[2][1] * mat[3][2]
      - mat[0][3] * mat[1][1] * mat[2][2] * mat[3][0]
      - mat[0][3] * mat[1][2] * mat[2][0] * mat[3][1];

  // not invertible
  if (detrminant == 0.0) {
    return std::vector<std::vector<double> >();
  }
  // inverse
  determinant = 1.0 / determinant;
  for (int i = 0; i < inverse.size(); ++i) {
    for (int j = 0; j < inverse[i].size(); ++j) {
      inverse[i][j] *= determinant;
    }
  }
  

  return inverse;
}




/*******************************End of decimation functions*****************************************/




/**********************************************************************************************************/


/**************************************** control_cb() *******************/
/* GLUI control callback                                                 */

void control_cb( int control )
{


// std::cout<< "open_file: "<<OPEN_FILE<<"   "<<"save_file: "<<OUTPUT_FILE<<"\n";

// std::cout<< "load_file: "<<LOAD_MESH<<"   "<<"save: "<<SAVE_FILE<<"\n";
if(control == SHADDING_ID){
  switch (curr_string){

        //flat shaded
        case 0 : 
                glShadeModel(GL_FLAT);
                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
               
                break;
        //smooth shaded
        case 1 : 
                glShadeModel(GL_SMOOTH);
                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        
                break;
        //wireframe
        case 2: 
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

                break;
        //shaded with mesh
        case 3:
        std::cout<<"print 1\n";
                glShadeModel(GL_SMOOTH);
                glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
                

      }
    }

  if ( control == LIGHT0_ENABLED_ID ) {
    if ( light0_enabled ) {
      glEnable( GL_LIGHT0 );
      light0_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT0 ); 
      light0_spinner->disable();
    }
  }
  else if ( control == LIGHT1_ENABLED_ID ) {
    if ( light1_enabled ) {
      glEnable( GL_LIGHT1 );
      light1_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT1 ); 
      light1_spinner->disable();
    }
  }
  else if ( control == LIGHT0_INTENSITY_ID ) {
    float v[] = { 
      light0_diffuse[0],  light0_diffuse[1],
      light0_diffuse[2],  light0_diffuse[3] };
    
    v[0] *= light0_intensity;
    v[1] *= light0_intensity;
    v[2] *= light0_intensity;

    glLightfv(GL_LIGHT0, GL_DIFFUSE, v );
  }
  else if ( control == LIGHT1_INTENSITY_ID ) {
    float v[] = { 
      light1_diffuse[0],  light1_diffuse[1],
      light1_diffuse[2],  light1_diffuse[3] };
    
    v[0] *= light1_intensity;
    v[1] *= light1_intensity;
    v[2] *= light1_intensity;

    glLightfv(GL_LIGHT1, GL_DIFFUSE, v );
  }
  else if ( control == ENABLE_ID )
  {
    glui2->enable();
  }
  else if ( control == DISABLE_ID )
  {
    glui2->disable();
  }
  else if ( control == SHOW_ID )
  {
    glui2->show();
  }
  else if ( control == HIDE_ID )
  {
    glui2->hide();
  }


  else if (control == OPEN_FILE)
  {
    strcpy(open_filename,"./");
    strcat(open_filename,open_filetext);
    strcat(open_filename,".smf");

    std::cout<< "open file name: "<< open_filename<<"\n";
  }
  else if (control == LOAD_MESH)
  {
smf.loadFile(open_filename);
glutPostRedisplay();

  }
  else if (control == OUTPUT_FILE)
  {
  strcpy(save_filename,"./");
  strcat(save_filename,save_filetext);
  strcat(save_filename,".smf");

  // std::cout<< "save file name: "<< save_filename<<"\n";

  }

  else if (control == SAVE_FILE)
  {
smf.saveFile(save_filename);
  }

  else if (control == DECIMATE){
    std::cout << "Value of k (choosing randomly amongst k candidates): " << k<< " edges. Collapsing  " << collapse_edges << " edges" << std::endl;
    smf.decimate(atoi(k), atoi(collapse_edges));
    glutPostRedisplay();

  }

}

/**************************************** myGlutKeyboard() **********/

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
  case 27: 
  case 'q':
    exit(0);
    break;
  };
  
  glutPostRedisplay();
}


/***************************************** myGlutMenu() ***********/

void myGlutMenu( int value )
{
  myGlutKeyboard( value, 0, 0 );
}


/***************************************** myGlutIdle() ***********/

void myGlutIdle()
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
  if ( glutGetWindow() != main_window ) 
    glutSetWindow(main_window);  

  /*  GLUI_Master.sync_live_all();  -- not needed - nothing to sync in this
                                       application  */

  glutPostRedisplay();
}

/***************************************** myGlutMouse() **********/

void myGlutMouse(int button, int button_state, int x, int y )
{
}


/***************************************** myGlutMotion() **********/

void myGlutMotion(int x, int y )
{
  glutPostRedisplay(); 
}

/**************************************** myGlutReshape() *************/

void myGlutReshape( int x, int y )
{
  int tx, ty, tw, th;
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
  glViewport( tx, ty, tw, th );

  xy_aspect = (float)tw / (float)th;

  glutPostRedisplay();
}


/************************************************** draw_axes() **********/
/* Disables lighting, then draws RGB axes                                */

void draw_axes( float scale )
{
  glDisable( GL_LIGHTING );

  glPushMatrix();
  glScalef( scale, scale, scale );

  glBegin( GL_LINES );
 
  glColor3f( 1.0, 0.0, 0.0 );
  glVertex3f( .8f, 0.05f, 0.0 );  glVertex3f( 1.0, 0.25f, 0.0 ); /* Letter X */
  glVertex3f( 0.8f, .25f, 0.0 );  glVertex3f( 1.0, 0.05f, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 1.0, 0.0, 0.0 ); /* X axis      */

  glColor3f( 0.0, 1.0, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 1.0, 0.0 ); /* Y axis      */

  glColor3f( 0.0, 0.0, 1.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 0.0, 1.0 ); /* Z axis    */
  glEnd();

  glPopMatrix();

  glEnable( GL_LIGHTING );
}


/***************************************** myGlutDisplay() *****************/

void myGlutDisplay()
{
  glClearColor( .9f, .9f, .9f, 1.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum( -xy_aspect*.04, xy_aspect*.04, -.04, .04, .1, 15.0 );

  glMatrixMode( GL_MODELVIEW );

  glLoadIdentity();
  glMultMatrixf( lights_rotation );
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  
  glLoadIdentity();
  glTranslatef( 0.0, 0.0, -2.6f );
  glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
 

  glScalef( scale, scale, scale );

  /***   These are _live_ variables, which are transparently 
    updated by GLUI ***/

  glPushMatrix();
  glTranslatef( 0.0, -0.25, 0.0 );
  glMultMatrixf( mesh_rotate );

  if(show_mesh){

    // glutSolidTorus( .15,.3,16,segments );
   

     switch (curr_string){

        //flat shaded
        case 0 : 
                glPushMatrix();
                glTranslatef( -.5, 0.0, 0.0 );
                glMultMatrixf( mesh_rotate );
                glColor3f(0.9f, 0.9f, 0.9f);
                glShadeModel(GL_FLAT);
          glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
                smf.display();
                glPopMatrix();
                // std::cout<<"check-------flat shaded-------\n";
                
                break;
        //smooth shaded
        case 1 : 
            glPushMatrix();
                glTranslatef( -.5, 0.0, 0.0 );
                glMultMatrixf( mesh_rotate );
                glColor3f(0.9f, 0.9f, 0.9f);
                glShadeModel(GL_SMOOTH);
                glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
                smf.display();
                glPopMatrix();
                // std::cout<<"check-------smooth shaded-------\n";
                break;
        //wireframe
        case 2: 
          
    glPushMatrix();
    glTranslatef( -.5, 0.0, 0.0 );
    glMultMatrixf( mesh_rotate );
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.9f, 0.9f, 0.9f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    smf.display();
    // std::cout<<"check-------wireframe-------\n";
    glDisable(GL_COLOR_MATERIAL);
    glPopMatrix();
    break;
        //shaded with mesh
        case 3:
             
    glPushMatrix();
    glTranslatef( -.5, 0.0, 0.0 );
    glMultMatrixf( mesh_rotate );
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glColor3f(0.9f, 0.9f, 0.9f);
    smf.display();
    // std::cout<<"print 2\n";
    // std::cout<<"check-------mesh and  shaded-------\n";
    glPopMatrix();
    glPushMatrix();
    glTranslatef( -.5, 0.0, 0.0 );
    glMultMatrixf( mesh_rotate );
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(0.0f, 0.0f, 0.0f);
    smf.display();
    glColor3f(0.9f, 0.9f, 0.9f);
    glDisable(GL_COLOR_MATERIAL);
    glPopMatrix();
    
    break;


      }

  }
  if ( show_axes )
    draw_axes(.52f);
  glPopMatrix();

  if ( show_text ) 
  {
    glDisable( GL_LIGHTING );  /* Disable lighting while we render text */
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( 0.0, 100.0, 0.0, 100.0  );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glColor3ub( 0, 0, 0 );
    glRasterPos2i( 10, 10 );

    /*  printf( "text: %s\n", text );              */

    /*** Render the live character array 'text' ***/
    int i;
    for( i=0; i<(int)strlen( string_list[curr_string] ); i++ )
      glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, string_list[curr_string][i] );
  }

  glEnable( GL_LIGHTING );


  glutSwapBuffers(); 
}


/**************************************** main() ********************/

int main(int argc, char* argv[])
{
  /****************************************/
  /*   Initialize GLUT and create window  */
  /****************************************/
 

  std::cout << smf;
// ~Smf();

  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 800, 600 );
 
  main_window = glutCreateWindow( "SMF" );

  
  glutDisplayFunc( myGlutDisplay );

  GLUI_Master.set_glutReshapeFunc( myGlutReshape );  
  GLUI_Master.set_glutKeyboardFunc( myGlutKeyboard );
  GLUI_Master.set_glutSpecialFunc( NULL );
  GLUI_Master.set_glutMouseFunc( myGlutMouse );
  glutMotionFunc( myGlutMotion );

  
 
  /****************************************/
  /*       Set up OpenGL lights           */
  /****************************************/

  glEnable(GL_LIGHTING);
  glEnable( GL_NORMALIZE );

  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  /****************************************/
  /*          Enable z-buferring          */
  /****************************************/

  glEnable(GL_DEPTH_TEST);

  /****************************************/
  /*         Here's the GLUI code         */
  /****************************************/

  printf( "GLUI version: %3.2f\n", GLUI_Master.get_version() );

  /*** Create the side subwindow ***/
  glui = GLUI_Master.create_glui_subwindow( main_window, 
              GLUI_SUBWINDOW_RIGHT );

  // /*** Add another rollout ***/
  GLUI_Rollout *options = new GLUI_Rollout(glui, "Options", false );
  new GLUI_Checkbox( options, "Show mesh", &show_mesh );
  new GLUI_Checkbox( options, "Show axes", &show_axes );
  new GLUI_Checkbox( options, "Show text", &show_text );

  /**** Add listbox ****/
  new GLUI_StaticText( glui, "" );
  GLUI_Listbox *list = new GLUI_Listbox( glui, "Mesh Options", &curr_string,SHADDING_ID,control_cb );
  int i;
  for( i=0; i<4; i++ )
    list->add_item( i, string_list[i] );

  new GLUI_StaticText( glui, "" );


  new GLUI_EditText(glui, "Open File:", GLUI_EDITTEXT_TEXT, open_filetext,OPEN_FILE,control_cb);
  new GLUI_Button(glui, "Load", LOAD_MESH, control_cb);

  
  new GLUI_EditText(glui, "Output File:", GLUI_EDITTEXT_TEXT, save_filetext,OUTPUT_FILE,control_cb);
  new GLUI_Button(glui,"Save", SAVE_FILE,control_cb);

new GLUI_StaticText( glui, "" );
new GLUI_StaticText( glui, "" );

new GLUI_EditText(glui, "Value for k: ", GLUI_EDITTEXT_TEXT, k ,K_VALUE,control_cb);
new GLUI_EditText(glui, "count to be removed:", GLUI_EDITTEXT_TEXT, collapse_edges,C_EDGES,control_cb);

new GLUI_Button(glui, "decimate", DECIMATE, control_cb );

new GLUI_StaticText( glui, "" );
new GLUI_StaticText( glui, "" );



  /****** A 'quit' button *****/

  
  new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );




  /**** Link windows to GLUI, and register idle callback ******/
  
  glui->set_main_gfx_window( main_window );


  /*** Create the bottom subwindow ***/
  glui2 = GLUI_Master.create_glui_subwindow( main_window, 
                                             GLUI_SUBWINDOW_BOTTOM );
  glui2->set_main_gfx_window( main_window );


  GLUI_Rotation *tor_rot = new GLUI_Rotation(glui2, "Mesh", mesh_rotate );
  tor_rot->set_spin( .98 );
  new GLUI_Column( glui2, false );


  GLUI_Translation *trans_xy = 
    new GLUI_Translation(glui2, "Objects XY", GLUI_TRANSLATION_XY, obj_pos );
  trans_xy->set_speed( .005 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_x = 
    new GLUI_Translation(glui2, "Objects X", GLUI_TRANSLATION_X, obj_pos );
  trans_x->set_speed( .005 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_y = 
    new GLUI_Translation( glui2, "Objects Y", GLUI_TRANSLATION_Y, &obj_pos[1] );
  trans_y->set_speed( .005 );
  new GLUI_Column( glui2, false );
  GLUI_Translation *trans_z = 
    new GLUI_Translation( glui2, "Objects Z", GLUI_TRANSLATION_Z, &obj_pos[2] );
  trans_z->set_speed( .005 );

  new GLUI_Column( glui2, false );

  new GLUI_Column( glui2, false );
  


#if 0
  /**** We register the idle callback with GLUI, *not* with GLUT ****/
  GLUI_Master.set_glutIdleFunc( myGlutIdle );
#endif
  

  /**** Regular GLUT main loop ****/
  
  glutMainLoop();



  return EXIT_SUCCESS;
}

