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
class Re_numbering
{
public:
  Re_numbering ();
  void run ();
private:
  void make_grid ();
  void distribute_dofs ();
  void output_results () const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;
  std::vector<typename DoFHandler<2>::active_cell_iterator> global_cells;
};

Re_numbering::Re_numbering ()
  :
  fe (1),
  dof_handler (triangulation)
{}

void  make_grid ()
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
void Re_numbering::distribute_dofs ()
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
 for (unsigned int i=0; i<dofs_per_cell; ++i )
 std::cout << "dof_indices in global mesh BEFORE DoFRenumbering "<< local_dof_indices[i]  << std::endl;
}

//........................

std::vector<types::global_dof_index> renum_local_dof_indices (dofs_per_cell);

DoFHandler<2>::active_cell_iterator
  celll = dof_handler.begin_active(),
  endcl = dof_handler.end();
  for (; celll!=endcl; ++celll)
  global_cells.push_back(celll);

   DoFRenumbering::cell_wise(dof_handler,global_cells);

 for (unsigned int i=0; i<global_cells.size(); ++i)
  {
    global_cells[i]->get_dof_indices (renum_local_dof_indices);
    for (unsigned int j=0; j<dofs_per_cell; ++j ){
     std::cout << "dof_indices in global mesh AFTER DoFRenumbering"<< renum_local_dof_indices[j]  << std::endl ;
    }
  }
}
//.............................................................................................
void Re_numbering::output_results () const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler (dof_handler);
  // vector of local_dof_indices for each cell
  data_out.build_patches ();
}

//.............................................................................................
void Re_numbering::run ()
{
  make_grid ();
  distribute_dofs ();
  output_results ();
}
//.............................................................................................
int main ()
{
  Re_numbering cell_wise_renumber;
  cell_wise_renumber.run ();
  return 0;
}
