#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
///#include <Eigen/Core>
#include <stdlib.h>
#include<numeric>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "../../../class_CASM_289H/eigen-git-mirror/Eigen/Dense"
#include "../../../class_CASM_289H/eigen-git-mirror/Eigen/LU"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////
// Defines three lattice vectors of a crystal
class Lattice
{
private:
    Eigen::Matrix3d LAT;

public:
    // different intilization instructors for different situations
    Lattice(Eigen::Matrix3d& M) { LAT = M; }

    Lattice(const Eigen::Vector3d& V1, const Eigen::Vector3d& V2, const Eigen::Vector3d& V3)
    {
        LAT.col(0) = V1;
        LAT.col(1) = V2;
        LAT.col(2) = V3;
    }
    // property
    Eigen::Vector3d Lattice_vector(int i) { return LAT.col(i); }
    void update(Eigen::Matrix3d M) { LAT=M; }
    Eigen::Matrix3d array() { return LAT; }


};

////////////////////////////////////////////////////////////////////////////////////////////
// Defines Cartesian position in a crystal
class Coordinate
{

public:
    Coordinate(Eigen::Vector3d coord)
    {
        //	std::vector<double> my_coord(coord);
        my_coord = coord;
    }

    Eigen::Vector3d get_coordinate() { return my_coord; }

    double get_x() { return my_coord(0); }
    double get_y() { return my_coord(1); }
    double get_z() { return my_coord(2); }

private:
    Eigen::Vector3d my_coord; // does this have to be public if it's used for the constructor?
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines position and type of atom in a crystal
class Site
{
public:
    Site(std::string atom_name, Coordinate& init_coord) : my_coord(init_coord), atom(atom_name) {}
    std::string get_atom() { return atom; }

    Eigen::Vector3d get_coordinate() { return my_coord.get_coordinate(); }

    // std::map<int, std::vector<double>> get_site(int i, Coordinate my_coord)
    //{

    //	std::vector<double> coord=my_coord.get_coordinate();
    //	std::map<int, std::vector<double>> init_site;
    // iit->first; //not sure what to put here
    // it->second;
    //	init_site.insert(std::make_pair(i, coord));
    //	return init_site;
    //}

private:
    std::string atom;
    Coordinate my_coord;
};

/////////////////////////////////////////////////////////////////////////////////////////
Site bring_within(const Site& site_outside_of_unit, const Lattice& unit_cell)
{
};
/////////////////////////////////////////////////////////////////////////////////////////
//Defines Cartesian matrix and translation vector
//of symmetry operation
class SymOp
{
};
/////////////////////////////////////////////////////////////////////////////////////////
//Defines Lattice and basis (collection of Site)
//in a crystal
class Structure
{
private:
    Lattice LatticeV;
    std::vector<Site> SitesV;

public:
    // different intilization instructors for different situations
    Structure(Lattice& StLattice, std::vector<Site>& StSites): LatticeV(StLattice), SitesV(StSites){}

};

///////////////////////////////////////////////////////////////////////////////////////
Structure read_poscar(const std::string& poscar_path)
{
	std::ifstream file(poscar_path);
    std::string title;
    std::getline(file,title);

    double scaling;
    file>>scaling;

    Eigen::Matrix3d lat_row_matrix;
    for(int i=0;i<3;i++)
    	for(int j=0;j<3;j++)
    		file>>lat_row_matrix(i,j);

    std::string species_line;
    //why needing two getline? because using >> to read data from file doesn't go to the next row!!!!but getline will go!
    std::getline(file,species_line);
    std::getline(file,species_line);
    std::istringstream iss(species_line);
	//defining iterator to read spaced words from file
    std::vector<std::string> species((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

    std::getline(file,species_line);
    std::istringstream iss2(species_line);
    //defining iterator to read spaced words from file
    std::vector<std::string> NumSpecies((std::istream_iterator<std::string>(iss2)), std::istream_iterator<std::string>());

    std::string coord_type;
    std::getline(file, coord_type);

    Eigen::Vector3d coord;
    std::vector<Eigen::Vector3d> raw_coordinate_values;
    while(!file.eof())
    {
    	for (int i=0;i<3;i++)
    		file>>coord(i);
    	raw_coordinate_values.push_back(coord);
    }

	//make Lattice
    Lattice latt(lat_row_matrix);
    //making the sites, creating one site per atom
    std::vector<Site> Sites;
    int t=0;
    for (int i=0;i<species.size();i++)
    	for (int j=0;j<std::stoi(NumSpecies[i]);j++)
    	{
    		Coordinate temp(raw_coordinate_values[t]);
    		Site Single_site(species[i],temp);
    		Sites.push_back(Single_site);
    		t=t+1;
    	}

    //making the structure
    Structure structure1(latt,Sites);
//for instance printing first row
    cout<<raw_coordinate_values[0][0]<<"\t"<<raw_coordinate_values[0][1]<<"\t"<<raw_coordinate_values[0][2]<<endl;
    return structure1;

};


// Defines a collection of Sites in a crystal
class Cluster
{
public:
    Cluster(std::vector<Site> sites) : my_sites(sites) {}

    int cluster_size() { return my_sites.size(); }

    Site get_site(int i) { return my_sites.at(i); }

    // std::map<int, std::vector<Site>> get_cluster(int i, std::vector<Site> my_sites)
    //{
    //     std::map<int, std::vector<Site>> my_cluster;
    // std::map<int, std::vector<Site>>::itterator it=my_cluster.begin();
    // for (int i=0; i<my_sites.size(); i++)
    //{
    //   my_cluster.insert(std::pair<int, std::vector<Site>>(i, my_sites));
    //}
    // return my_cluster;
    //}

private:
    std::vector<Site> my_sites;
};
/*
struct SiteCompare_f
{
    bool operator()(const Site& other)
    {
    }
};
struct ClusterCompare_f
{
    bool operator()(const Cluster& other)
    {
    }
};
Site operator*(const Site& site, const SymOp& transformation)
{
}
Cluster operator*(const Cluster& site, const SymOp& transformation)
{
}
std::vector<SymOp> make_factor_group(const Structure& struc)
{
}
*/
int main()
{
//testing of Lattice


Eigen::Matrix3d O=Eigen::Matrix3d::Ones(3,3);

Lattice M2(O);
cout<<M2.Lattice_vector(1)<<endl;

//testing of read_poscar


//read_poscar("POSCAR.txt");


read_poscar("POSCAR");





return 0;

}
