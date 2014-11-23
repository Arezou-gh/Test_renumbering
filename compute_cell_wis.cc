// To test how cell_wise renumbering works

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

using namespace dealii;
class compute_cell_wis
{
public:
  compute_cell_wis ();
  void run ();
private:
  void make_grid ();
  void distribute_dofs ();
  void output_results () const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
   DoFHandler<2>        dof_handler;
  const std::vector<typename DoFHandler<2>::active_cell_iterator> global_cells;
};

compute_cell_wis::compute_cell_wis ()
  :
  fe (1),
  dof_handler (triangulation)
{}

void compute_cell_wis::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  std::ofstream out ("grid.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
  std::cout << "Grid written to grid.eps" << std::endl;

   std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
}
//.............................................................................................
void compute_cell_wis::distribute_dofs ()
{
 dof_handler.distribute_dofs (fe);
 std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
const unsigned int   dofs_per_cell = fe.dofs_per_cell;

//......................

std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

 DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
{
cell->get_dof_indices (local_dof_indices);
std::cout << "dof_indices in global mesh BEFORE DoFRenumbering ";
std::cout<<cell->center()(0)<<" "<<cell->center()(1)<<":";
 for (unsigned int i=0; i<dofs_per_cell; ++i )
  std::cout<<" "<< local_dof_indices[i];
  std::cout<< std::endl;
}
std::vector< types::global_dof_index >  renumbering;
std::vector< types::global_dof_index >  inverse_renumbering;

DoFRenumbering::compute_cell_wis( renumbering, inverse_renumbering, dof_handler,global_cells);

for (unsigned int j=0; j<renumbering.size() ; ++j){
std::cout<< " renumbering vector :"  <<renumbering[j];
std::cout  << std::endl ;
}

for (unsigned int k=0; k<inverse_renumbering.size() ; ++k){
std::cout<< " inverse_renumbering vector: " <<inverse_renumbering[k];
std::cout  << std::endl ;
}

//........................
/*
std::vector<types::global_dof_index> renum_local_dof_indices (dofs_per_cell);

 std::vector<typename DoFHandler<2>::active_cell_iterator> tmp_cells;

DoFHandler<2>::active_cell_iterator
  celll = dof_handler.begin_active(),
  endcl = dof_handler.end();
  for (; celll!=endcl; ++celll)
  tmp_cells.push_back(celll);

global_cells.resize(tmp_cells.size());
for (unsigned int i=0; i<tmp_cells.size(); ++i)
global_cells[tmp_cells.size()-1-i] = tmp_cells[i];

   DoFRenumbering::cell_wise(dof_handler,global_cells);
cell = dof_handler.begin_active(),
  endc = dof_handler.end();

for (; cell!=endc; ++cell)
  {
    cell->get_dof_indices (renum_local_dof_indices);
 std::cout << "dof_indices in global mesh AFTER DoFRenumbering"; 
std::cout<<cell->center()(0)<<" "<<cell->center()(1)<<":";  
for (unsigned int j=0; j<dofs_per_cell; ++j )
     std::cout<<" " << renum_local_dof_indices[j];
std::cout  << std::endl ;
    
  }
*/
}
//.............................................................................................
void compute_cell_wis::output_results () const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  // vector of local_dof_indices for each cell
  data_out.build_patches ();
}

//.............................................................................................
void compute_cell_wis::run ()
{
  make_grid ();
  distribute_dofs ();
  output_results ();
}
//.............................................................................................
int main ()
{
  compute_cell_wis cell_wise_renumber;
  cell_wise_renumber.run ();
  return 0;
}
