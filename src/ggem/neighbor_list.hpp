#ifndef LOCAL_CELLS_HPP
#define LOCAL_CELLS_HPP

#include "mesh.h"
#include <map>
#include <vector>

//WEST  -x <-------------------> +x EAST
//SOUTH -y <-------------------> +y NORTH
//BACK  -z <-------------------> +z FRONT
enum class BoundaryName {EAST, WEST, NORTH, SOUTH, FRONT, BACK};

class NeighborList{
    public:
        NeighborList(){};

        NeighborList(double Lx, double Ly, double Lz, double rcut);

        ~NeighborList();

        void set_pbc(const int p[3]);

        // Recalculates all neighbors
        void build_list(const Mesh* msh);

        // Returns a pointer to a list of cells for vertex ivert
        int* get_cells_nbring_vert(const int ivert, int* nc);

        // Returns a pointer to a list of cells for boundary vertex
        // indexed by bvert_indx. bndry = 't' or 'b', represents top
        // or bottom boundary
        int* get_cells_nbring_bvert(const int bvert_indx[2],
                const char bndry, int* nc);

    private:
        //Sorts vertices into boxes
        void sort_verts(const int num_verts, const double* coordinates);

        void build_verts_nbring_vert_list();

        void NeighborList::build_cells_nbring_vert_list();

        //Returns a vector whose each element is a vector containing the
        //indices of the neighbors of a box
        void get_boxes_nbring_box(const std::vector<int> indx,
                std::vector<std::vector<int> >& nbrs);

        //Returns a vector whose each element is a vector containing the
        //indices of the boundary vertices incident to a boundary box
        void get_bverts_nbring_bbox(const std::vector<int>& box_indx,
                const char boundary, std::vector<std::vector<int> >& nbrs);

        bool is_nbring_point(const double* p1, const double* p2);

        //Returns true if a box_indx represents a boundary box
        bool is_boundary_box(const std::vector<int> box_indx, char* boundary);

        unsigned long int box_indx_to_ibox(const std::vector<int>& box_indx);

        void ibox_to_box_indx(const unsigned long int ibox,
                std::vector<int>& box_indx);

        unsigned long int bvert_indx_to_ibvert(const std::vector<int>& bvert_indx);

        void ibvert_to_bvert_indx(const unsigned long int ibvert,
                std::vector<int>& bvert_indx);

        //Applies periodic boundary condition to a box_indx or bvert_indx
        void apply_pbc_indx(const char indx_of, std::vector<int>& indx);


        double Lx_;
        double Ly_;
        double Lz_;
        double rcut_;
        int pbc_[3];

        int num_boxes_x_;
        int num_boxes_y_;
        int num_boxes_z_;

        std::vector<double> offset_;

        std::map<unsigned long int, std::vector<int> > box_to_verts_;

        std::vector<std::vector<int> > vert_to_verts_;

        std::vector<std::vector<int> > vert_to_cells_;

        //top boundary
        std::map<unsigned long int, std::vector<int> > tbvert_to_cells_;

        //bottom boundary
        std::map<unsigned long int, std::vector<int> > bbvert_to_cells_;
};

#endif // LOCAL_CELLS_HPP
