#include "local_cells.hpp"
#include <cmath>

/*****************************************************************************/

NeighborList::NeighborList(double Lx, double Ly, double Lz, double rcut)
    : Lx_(Lx), Ly_(Ly), Lz_(Lz), rcut_(rcut), pbc_ {0, 0, 0}
{
    num_boxes_x_ = lround(Lx_/rcut_) + 2; //Adding a ghost box on either side
    num_boxes_y_ = lround(Ly_/rcut_) + 2;
    num_boxes_z_ = lround(Lz_/rcut_) + 2;

    offset_ = std::vector<double>(3, rcut_);
}

/*****************************************************************************/

void NeighborList::set_pbc(const int p[3])
{
    pbc_[0] = p[0]; pbc_[1] = p[2]; pbc_[2];
}

/*****************************************************************************/

void NeighborList::build_list(Mesh* msh)
{
    //Sort vertices into boxes
    int num_verts = mesh_get_num_vertices(msh);
    double* coordinates = mesh_get_coordinates(msh, NULL);
    sort_verts(num_verts, coordinates);

    //Build list of vertices neighboring each vertex
    build_verts_nbring_vert_list();

    //Loop over all vertices and add cells
    for (auto& each : vert_to_cells_){
        each.clear(); //capacity may not change
    }

    int* conns;
    int nc;

    for (int ivert=0; ivert < num_verts; ++ivert){
        for (const int jvert : vert_to_verts_[ivert]){
            conns = mesh_get_conns(msh, 0, 2, jvert, &nc);
            for (int i=0; i < nc; ++i){
                vert_to_cells_[ivert].push_back(conns[i]);
            }
        }
    }

    //Remove duplicate cells for each vertex & release unused memory
    for (auto& each : vert_to_cells_){
        std::sort(each.begin(), each.end());
        auto last = std::unique(each.begin(), each.end());
        each.erase(last, each.end());
        each.shrink_to_fit();
    }

    //
    //Loop over all occupied boxes to check if any boundary box is occupied
    //Add additional logic to account for densely filled domain

    tbverts_to_cells.clear();
    bbverts_to_cells.clear();

    std::vector<int> bvert_indx(2,0);
    char boundary;
    bool is_bndry;

    for (auto& keyval : box_to_verts_){
        ibox = keyval.first;
        ibox_to_box_indx(ibox, box_indx);
        is_bndry = is_boundary_box(box_indx, &boundary);
        if (is_bndry){
            get_bverts_nbring_box(box_indx, boundary, nbrs);
            for (int& ivert : keyval.second){
                coords_ivert = mesh_get_coords(msh, ivert, NULL);
                for (){
                    coords_jvert = ;
                    if(is_nbring_point(coords_ivert, coords_jvert)){
                        conns = mesh_get_conns(msh,0, 2, jvert, &nc);
                        for (int i=0; i < nc; ++i){
                            tbvert_to_cells_[ibvert].push_back(conns[i]);
                        }
                    }
                }
            }
        }
    }

    //Remove duplicate cells for each boundary vertex & release unused memory
    std::vector<int> cells;
    for (auto& each : tbvert_to_cells_){
        cells = each.second;
        std::sort(cells.begin(), cells.end());
        auto last = std::unique(cells.begin(), cells.end());
        cells.erase(last, cells.end());
        cells.shrink_to_fit();
    }

}

/*****************************************************************************/

//Sort vertices into boxes
void NeighborList::sort_verts(const int num_verts,
        const double* coordinates);
{
    int ibox;
    std::vector<int> box_indx(3,0);
    double* coords_ivert;

    box_to_verts_.clear();

    for (int ivert=0; ivert < num_verts; ++ivert){
        coords_ivert = &coordinates[3*ivert];

        for (int i=0; i < 3; ++i){
            box_indx[i] = 1 + int(coords_ivert[i]/rcut_);
        }
        ibox = box_indx_to_ibox(box_indx);
        box_to_verts_[ibox].push_back(ivert);
    }
}

/*****************************************************************************/

void NeighborList::build_verts_nbring_vert_list()
{
    for (auto& i: vert_to_verts_){
        i.clear(); //capacity may not change 
    }

    //Loop over all boxes to find vertices within cutoff
    int ibox;
    int ibox_nbr;
    std::vector<int> box_indx(3,0);
    double* coords_ivert;
    double* coords_jvert;
    std::vector<std::vector<int> > box_nbrs(13, std::vector<int>(3));

    for (auto& keyval : box_to_verts_){
        ibox = keyval.first;
        ibox_to_box_indx(ibox, box_indx);
        get_box_nbrs(box_indx, box_nbrs);

        //For each vertex within a box add all cells incident to it
        for (const int& ivert : keyval.second){
            coords_ivert = mesh_get_coords(msh, ivert, NULL);
            for (const int& jvert : keyval.second){
                vert_to_verts_[ivert].push_back(jvert);
            }
            //Search in neighboring boxes
            for (auto& box_nbr : box_nbrs){
                ibox_nbr = box_indx_to_ibox(box_nbr);
                //If ibox_nbr is not empty
                if (box_to_verts_.count(ibox_nbr) == 1){
                    for (const int& jvert : box_to_verts_[ibox_nbr]){
                        coords_jvert = mesh_get_coords(msh, jvert, NULL);
                        if (is_nbring_point(coords_ivert, coords_jvert)){
                            vert_to_verts_[ivert].push_back(jvert);
                            vert_to_verts_[jvert].push_back(ivert);
                        }
                    }
                }
            }
        }
    }

    //Release unused memory
    for (auto& each: vert_to_verts_){
        each.shrink_to_fit();
    }
}

/*****************************************************************************/

void NeighborList::build_cells_nbring_vert_list()
{

}

/*****************************************************************************/

//Returns neighbor boxes of a given box
void NeighborList::get_boxes_nbring_box(const std::vector<int> indx,
        std::vector<std::vector<int> >& nbrs)
{
    int ix = indx[0];
    int iy = indx[1];
    int iz = indx[2];

    // Store Ist half
    // Sharing face (6)
    nbrs[0][0] = ix ;  nbrs[0][1] = iy  ;  nbrs[0][2] = iz+1;
    nbrs[1][0] = ix ;  nbrs[1][1] = iy+1;  nbrs[1][2] = iz  ;
    nbrs[2][0] = ix+1; nbrs[2][1] = iy  ;  nbrs[2][2] = iz  ;

    // Sharing a side(12)
    nbrs[3][0] = ix  ;  nbrs[3][1] = iy+1;  nbrs[3][2] = iz+1; 
    nbrs[4][0] = ix+1;  nbrs[4][1] = iy  ;  nbrs[4][2] = iz+1;
    nbrs[5][0] = ix-1;  nbrs[5][1] = iy  ;  nbrs[5][2] = iz+1;
    nbrs[6][0] = ix-1;  nbrs[6][1] = iy+1;  nbrs[6][2] = iz  ;
    nbrs[7][0] = ix+1;  nbrs[7][1] = iy+1;  nbrs[7][2] = iz  ;
    nbrs[8][0] = ix  ;  nbrs[8][1] = iy+1;  nbrs[8][2] = iz-1;

    // Sharing a corner (8)
    nbrs[ 9][0] = ix+1;  nbrs[ 9][0] = iy+1; nbrs[ 9][0] = iz+1;
    nbrs[10][0] = ix-1;  nbrs[10][0] = iy+1; nbrs[10][0] = iz+1;
    nbrs[11][0] = ix+1;  nbrs[11][0] = iy+1; nbrs[11][0] = iz-1;
    nbrs[12][0] = ix-1;  nbrs[12][0] = iy+1; nbrs[12][0] = iz-1;

    //Check for periodicity
    for (auto& i : nbrs){
        apply_pbc_box_indx(i);
    }
}

/*****************************************************************************/

//Returns neighbor boxes of a given vertex in boundary (top or bottom)
void NeighborList::get_bverts_nbring_bbox(const std::vector<int>& box_indx,
        const char bndry, std::vector<std::vector<int> >& nbrs)
{
    int ix = box_indx[0];
    int iz = box_indx[2];
    int iy = box_indx[1];

    nbrs[ 0][0] = ix-2; nbrs[ 0][2] = iz-1;
    nbrs[ 1][0] = ix-1; nbrs[ 1][2] = iz-1;
    nbrs[ 2][0] = ix;   nbrs[ 2][2] = iz-1;
    nbrs[ 3][0] = ix+1; nbrs[ 3][2] = iz-1;
    nbrs[ 4][0] = ix-2; nbrs[ 4][2] = iz;
    nbrs[ 5][0] = ix-1; nbrs[ 5][2] = iz;
    nbrs[ 6][0] = ix;   nbrs[ 6][2] = iz;
    nbrs[ 7][0] = ix+1; nbrs[ 7][2] = iz;
    nbrs[ 8][0] = ix-1; nbrs[ 8][2] = iz-2;
    nbrs[ 9][0] = ix-1; nbrs[ 9][2] = iz+1;
    nbrs[10][0] = ix;   nbrs[10][2] = iz-2;
    nbrs[11][0] = ix;   nbrs[11][2] = iz+1;

    if (c=='t'){
        for (int i=0; i < 12; ++i){
            nbrs[i][1] = num_boxes_y_;
        }
    else if (c=='b'){
        for (int i=0; i < 12; ++i){
            nbrs[i][1] = 1;
        }
    }

    //Check for periodicity
    for (auto& each : nbrs){
        apply_pbc_bvert_indx(i);
    }
}

/*****************************************************************************/

//Returns true if two points are within cutoff distance
bool NeighborList::is_nbring_point(const double* p1, const double* p2)
{
    double rcut_sq = rcut_*rcut_;

    double rxij = p2[0] - p1[0];
    double ryij = p2[1] - p1[1];
    double rzij = p2[2] - p1[2];

    //Minimum image convention, no-op if not PBC
    rxij -= round(rxij/Lx_)*Lx_;
    ryij -= round(ryij/Lx_)*Ly_;
    rzij -= round(rzij/Lx_)*Lz_;

    double rij_sq = rxij*rxij + ryij*ryij + rzij*rzij;

    if (rij_sq < rcut_sq) {
        return true;
    }
    else {
        return false;
    }
}

/*****************************************************************************/

//Returns true if a box_indx represents a boundary box
void NeighborList::is_boundary_box(std::vector<int>& box_indx, char* boundary)
{

}

/*****************************************************************************/

unsigned long int NeighborList::box_indx_to_ibox(const std::vector<int>& box_indx)
{

}

/*****************************************************************************/

void NeighborList::ibox_to_box_indx(const unsigned long int ibox,
        std::vector<int>& box_indx)
{

}

/*****************************************************************************/

unsigned long int NeighborList::bvert_indx_to_ibvert(
        const std::vector<int>& bvert_indx)
{

}

/*****************************************************************************/

void NeighborList::ibvert_to_bvert_indx(const unsigned long int ibvert,
                std::vector<int>& bvert_indx)
{

}

/*****************************************************************************/

//Applies periodic boundary condition to a box_indx or bvert_indx
void NeighborList::apply_pbc_indx(const char indx_of, std::vector<int>& indx)
{
    if (pbc_[0] == 1){
        if (indx[0] < 1){
            indx[0] += (num_boxes_x_-2);
        }
        else if {
            if (indx_of == 'b') {
                if (indx[0] > (num_boxes_x_-2)){
                    indx[0] -= (num_boxes_x_-2);
                }
            }
            else if (indx_of == 'v'){
                if (indx[0] > (num_boxes_x_-1)){
                    indx[0] -= (num_boxes_x_-2);
                }
            }
        }
    }

    if (pbc_[1] == 1){
        if (indx[1] < 1){
            indx[1] += (num_boxes_y_-2);
        }
        else if {
            if (indx_of == 'b'){
                if (indx[0] > (num_boxes_y_-2)){
                    indx[0] -= (num_boxes_y_-2);
                }
            }
            else if (indx_of == 'v'){
                if (indx[0] > (num_boxes_y_-1)){
                    indx[0] -= (num_boxes_y_-2);
                }
            }
        }
    }

    if (pbc_[2] == 1){
        if (indx[2] < 1){
            indx[2] += (num_boxes_z_-2);
        }
        else if {
            if (indx_of == 'b'){
                if (indx[0] > (num_boxes_z_-2)){
                    indx[0] -= (num_boxes_z_-2);
                }
            }
            else if (indx_of == 'v'){
                if (indx[0] > (num_boxes_z_-1)){
                    indx[0] -= (num_boxes_z_-2);
                }
            }
        }
    }
}

/*****************************************************************************/

int* NeighborList::get_cells_nbring_vert(const int ivert, int* nc)
{
    *nc = vert_to_cells_[ivert].size();
    return vert_to_cells_[ivert].data();
}

/*****************************************************************************/

//indx[0] = nx, indx[1] = nz
int* NeighborList::get_cells_nbring_bvert(const int indx[2],
        const char bndry, int* nc)
{
    *nc = 0;
    ret = NULL;

    if (bndry == 't'){
        std::vector<int> bvert_indx{indx[0], indx[1]};
        if (tbvert_to_cells_.count(ibvert)==1){
            *nc = tbvert_to_cells_[ibvert].size();
            ret  = tbvert_to_cells_[ibvert].data();
        }
    }
    else if (bndry == 'b'){
        if (bbvert_to_cells_.count(ibvert)==1){
            *nc = bbvert_to_cells_[indx].size();
            ret = bbvert_to_cells_[indx].data();
        }
    }
    return ret;
}

/*****************************************************************************/
