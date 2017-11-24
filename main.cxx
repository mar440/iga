#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
#include <vtkPolyVertex.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vector>
#include <math.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vector>
#include <vtkFieldData.h>



double length[] = {1.0, 1.0, 1.0};


int nEl_IGA[] = {10, 20, 5};
//int nEl_IGA_x = 10;
//int nEl_IGA_y = 20;
//int nEl_IGA_z = 5;


bool asciiOrBinaryVtu = true;
bool display_mesh = false;


using namespace std;

int main(int argc, char *argv[])
{

   // - geometry setting ()

    int cnt;
    int nEl = 1;

    std::string filename = "iga.vtu";//argv[1];

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

// -----------------------------------------------------------------------------------------


        int nKnot_x = nEl_IGA[0] + 1;
        int nKnot_y = nEl_IGA[1] + 1;
        int nKnot_z = nEl_IGA[2] + 1;
        int nCP_x = nKnot_x + 1;
        int nCP_y = nKnot_y + 1;
        int nCP_z = nKnot_z + 1;
        int nx = (6 + nEl_IGA[0] - 1);
        int ny = (6 + nEl_IGA[1] - 1);
        int nz = (6 + nEl_IGA[2] - 1);
        vector < double > t0(nx + ny + nz,1);


        t0[0          ] = 0;
        t0[1          ] = 0;
        t0[2          ] = 0;

        t0[0 + nx     ] = 0;
        t0[1 + nx     ] = 0;
        t0[2 + nx     ] = 0;

        t0[0 + nx + ny] = 0;
        t0[1 + nx + ny] = 0;
        t0[2 + nx + ny] = 0;

        double knot_step_x = 1.0 / double (nEl_IGA[0]);
        double knot_step_y = 1.0 / double (nEl_IGA[1]);
        double knot_step_z = 1.0 / double (nEl_IGA[2]);
        for (int i = 0; i < nEl_IGA[0] - 1; i++){
            t0[3           + i] = knot_step_x * (i + 1);
        }
        for (int i = 0; i < nEl_IGA[1] - 1; i++){
            t0[3 + nx      + i] = knot_step_y * (i + 1);
        }
        for (int i = 0; i < nEl_IGA[2] - 1; i++){
            t0[3 + nx + ny + i] = knot_step_z * (i + 1);
        }


        for (int i = 0; i < nx; i++)
            cout << t0[i] << " ";
        cout << endl;
        for (int i = 0; i < ny; i++)
            cout << t0[i + nx] << " ";
        cout << endl;
        for (int i = 0; i < nz; i++)
            cout << t0[i + nx + ny] << " ";
        cout << endl;



    vtkSmartPointer<vtkDoubleArray> _vtkDataArray0 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray0->SetName("Knots_0");
    _vtkDataArray0->SetNumberOfComponents(1);
    _vtkDataArray0->SetNumberOfTuples(t0.size());
    for (int i = 0; i < t0.size(); i++){
        _vtkDataArray0->SetTuple(i, t0.data() + i);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray0);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkDoubleArray> _vtkDataArray1 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray1->SetName("Weights_0");
    _vtkDataArray1->SetNumberOfComponents(1);
    vector < double > t1(nCP_x * nCP_y * nCP_z,1);
    _vtkDataArray1->SetNumberOfTuples(t1.size());


    for (int i = 0; i < t1.size(); i++){
        _vtkDataArray1->SetTuple(i, t1.data() + i);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray1);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray2 = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray2->SetName("PatchData_0");
    _vtkDataArray2->SetNumberOfComponents(1);
    vector < int > t2(  { 2,  2,  2, nCP_x,  nCP_y,  nCP_z,
            6 + nEl_IGA[0] - 1,  6 + nEl_IGA[1] - 1,   6 + nEl_IGA[2] - 1 } );

    for (int i = 0; i < t2.size(); i++){
        _vtkDataArray2->InsertValue(i,t2[i]);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray2);
// -----------------------------------------------------------------------------------------

    int nP_IGA = t2[3] * t2[4] * t2[5];
// number of IGA nodes
    vtkSmartPointer<vtkIntArray> _vtkDataArray3 = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray3->SetName("RegionId");
    _vtkDataArray3->SetNumberOfComponents(1);
    vector < int > t3();
    _vtkDataArray3->SetNumberOfTuples(nP_IGA);

    for (int i = 0; i < nP_IGA; i++){
        _vtkDataArray3->InsertValue(i, 0);
    }
//    for (int i = nP_IGA; i < nP_IGA; i++){
//        _vtkDataArray3->InsertValue(i, 1);
//    }
    unstructuredGrid->GetPointData()->AddArray(_vtkDataArray3);
// -----------------------------------------------------------------------------------------

    double _x,_y,_z;
    int nPoints = 0;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vector < double > vPointsIGA(3 * nP_IGA);
    double dxyzIGA[3];
    double shiftIGA[3];
    for (int i = 0 ; i < 3 ; i++){
        dxyzIGA[i] = length[0] / double (nEl_IGA[i] + 1) ;
        shiftIGA[i] = - 0.5 * length[i];
    }

    cnt = 0;
    for (int kk = 0; kk <  nCP_z; kk++){
        for (int jj = 0; jj <  nCP_y; jj++){
            for (int ii = 0; ii <  nCP_x; ii++){
                _x = ii * dxyzIGA[0] + shiftIGA[0];
                _y = jj * dxyzIGA[1] + shiftIGA[1];
                _z = kk * dxyzIGA[2] + shiftIGA[2];
                vPointsIGA[cnt + 0] = _x;
                vPointsIGA[cnt + 1] = _y;
                vPointsIGA[cnt + 2] = _z;
                cnt += 3;
            }
        }
    }


    for (int i = 0; i <  nP_IGA ; i++){
        _x = vPointsIGA[3 * i + 0];
        _y = vPointsIGA[3 * i + 1];;
        _z = vPointsIGA[3 * i + 2];;
        points->InsertNextPoint(_x, _y, _z);
        nPoints++;
    }



    vtkSmartPointer<vtkIdList> plvx_ids = vtkSmartPointer<vtkIdList>::New();



    vector < int > v_plVert( nP_IGA );
    for (int i = 0 ; i < nP_IGA ; i++)
        v_plVert[i] = i;

    for (int i = 0 ; i < v_plVert.size(); i++){
        plvx_ids->InsertId(i,v_plVert[i]);
    }
    unstructuredGrid->Allocate(nEl,1);
    unstructuredGrid->InsertNextCell (VTK_POLY_VERTEX, plvx_ids);
    unstructuredGrid->SetPoints(points);



    int int_PartId = 0;
    int int_MatId = 1;
    int int_FormId = 900;
    int int_PieceId = 1;
    int int_RegionId = 0;

// PartitionId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PartId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PartId->SetName("PartitionId");
    _vtkDataArray_PartId->SetNumberOfComponents(1);
    _vtkDataArray_PartId->InsertValue(0,int_PartId);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PartId);
// MaterialId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_MatId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_MatId->SetName("MaterialId");
    _vtkDataArray_MatId->SetNumberOfComponents(1);
    _vtkDataArray_MatId->InsertValue(0,int_MatId);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_MatId);
// FormulationId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_FormId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_FormId->SetName("FormulationId");
    _vtkDataArray_FormId->SetNumberOfComponents(1);
    _vtkDataArray_FormId->InsertValue(0, int_FormId);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_FormId);

// PieceId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PieceId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PieceId->SetName("PieceId");
    _vtkDataArray_PieceId->SetNumberOfComponents(1);
    _vtkDataArray_PieceId->InsertValue(0,int_PieceId);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PieceId);
// RegionId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_RegionId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_RegionId->SetName("RegionId");
    _vtkDataArray_RegionId->SetNumberOfComponents(1);
    _vtkDataArray_RegionId->InsertValue(0, int_RegionId);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_RegionId);

#if 1

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());

    if (asciiOrBinaryVtu){
      writer->SetDataModeToAscii();
    }
    else{
      writer->SetDataModeToBinary();
    }


#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  writer->Write();


  cout << "number of Points:     " << nPoints << endl;


  // Read and display file for verification that it was written correclty
    if (display_mesh){
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
          vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(filename.c_str());
        reader->Update();

        vtkSmartPointer<vtkDataSetMapper> mapper =
          vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());

        vtkSmartPointer<vtkActor> actor =
          vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        vtkSmartPointer<vtkRenderer> renderer =
          vtkSmartPointer<vtkRenderer>::New();
        vtkSmartPointer<vtkRenderWindow> renderWindow =
          vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
          vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        renderer->AddActor(actor);
        renderer->SetBackground(.3, .6, .3); // Background color green

        renderWindow->Render();
        renderWindowInteractor->Start();
    }

#endif
  return EXIT_SUCCESS;
}

