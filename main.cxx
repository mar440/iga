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
#include <vtkCellData.h>
#include <vector>
#include <vtkFieldData.h>

int nEx = 6;
int nEy = 3;
int nEz = 4;
int nSx = 1;
int nSy = 2;
int nSz = 1;

int nEl_IGA_x = 5;


bool asciiOrBinaryVtu = true;
bool display_mesh = false;


using namespace std;

int main(int argc, char *argv[])
{

   // - geometry setting ()
    double length[] = {4.0, 2.0, 2.0};

    int cnt;
    int nEl = 1;

    std::string filename = "iga.vtu";//argv[1];

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkDoubleArray> _vtkDataArray0 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray0->SetName("Knots_0");
    _vtkDataArray0->SetNumberOfComponents(1);


        int nKnot_x = nEl_IGA_x + 1;
        int nCP_x = nKnot_x + 1;
        int n0 = (6 + nEl_IGA_x - 1);
        vector < double > t0( 3 * n0,1);


        t0[0 + n0 * 0] = 0;
        t0[1 + n0 * 0] = 0;
        t0[2 + n0 * 0] = 0;

        t0[0 + n0 * 1] = 0;
        t0[1 + n0 * 1] = 0;
        t0[2 + n0 * 1] = 0;

        t0[0 + n0 * 2] = 0;
        t0[1 + n0 * 2] = 0;
        t0[2 + n0 * 2] = 0;

        double knot_step = 1.0 / double (nEl_IGA_x);
        for (int i = 0; i < nEl_IGA_x - 1; i++){
            t0[3 + n0 * 0 + i] = knot_step * (i + 1);
            t0[3 + n0 * 1 + i] = knot_step * (i + 1);
            t0[3 + n0 * 2 + i] = knot_step * (i + 1);
        }


    for (int j = 0; j < 3; j++){
        for (int i = 0; i < n0; i++){
        cout << t0[j * n0 + i] << " ";
        }
        cout << endl;
    }



    _vtkDataArray0->SetNumberOfTuples(t0.size());
    double tuple[] = {0};
    for (int i = 0; i < t0.size(); i++){
        tuple[0] = t0[i];
        _vtkDataArray0->SetTuple(i, tuple);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray0);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkDoubleArray> _vtkDataArray1 = vtkSmartPointer<vtkDoubleArray>::New();
    _vtkDataArray1->SetName("Weights_0");
    _vtkDataArray1->SetNumberOfComponents(1);
    vector < double > t1(nCP_x * nCP_x * nCP_x,1);


    _vtkDataArray1->SetNumberOfTuples(t1.size());


    for (int i = 0; i < t1.size(); i++){
        tuple[0] = t1[i];
        _vtkDataArray1->SetTuple(i, tuple);
    }
    unstructuredGrid->GetFieldData()->AddArray(_vtkDataArray1);
// -----------------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray2 = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray2->SetName("PatchData_0");
    _vtkDataArray2->SetNumberOfComponents(1);
    vector < int > t2(  { 2,  2,  2,
                            nCP_x,  nCP_x,  nCP_x,
                            6 + nEl_IGA_x - 1,  6 + nEl_IGA_x - 1,   6 + nEl_IGA_x - 1 } );
    _vtkDataArray2->SetNumberOfTuples(t2.size());


    for (int i = 0; i < t2.size(); i++){
        tuple[0] = t2[i];
        _vtkDataArray2->SetTuple(i, tuple);
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
        tuple[0] = 0;
        _vtkDataArray3->SetTuple(i, tuple);
    }
    for (int i = nP_IGA; i < nP_IGA; i++){
        tuple[0] = 1;
        _vtkDataArray3->SetTuple(i, tuple);
    }
    unstructuredGrid->GetPointData()->AddArray(_vtkDataArray3);
// -----------------------------------------------------------------------------------------

    double _x,_y,_z;
    int nPoints = 0;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vector < double > vPointsIGA(3 * nP_IGA);
    double dxyzIGA[3];
    double shiftIGA[3];
    for (int i = 0 ; i < 3 ; i++){
        dxyzIGA[i] = length[i] / double (nEl_IGA_x + 1) ;
        shiftIGA[i] = - 0.5 * length[i];
    }
    shiftIGA[0] = -length[0];

    cnt = 0;
    for (int kk = 0; kk <  nCP_x; kk++){
        for (int jj = 0; jj <  nCP_x; jj++){
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
    vtkSmartPointer<vtkIdList> hexa_ids = vtkSmartPointer<vtkIdList>::New();



    vector < int > v_plVert( nP_IGA );
    for (int i = 0 ; i < nP_IGA ; i++)
        v_plVert[i] = i;

    for (int i = 0 ; i < v_plVert.size(); i++){
        plvx_ids->InsertId(i,v_plVert[i]);
    }
    unstructuredGrid->Allocate(nEl,1);
    unstructuredGrid->InsertNextCell (VTK_POLY_VERTEX, plvx_ids);
    unstructuredGrid->SetPoints(points);




// PartitionId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PartId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PartId->SetName("PartitionId");
    _vtkDataArray_PartId->SetNumberOfComponents(1);
    _vtkDataArray_PartId->SetNumberOfTuples(nEl);
// MaterialId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_MatId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_MatId->SetName("MaterialId");
    _vtkDataArray_MatId->SetNumberOfComponents(1);
    _vtkDataArray_MatId->SetNumberOfTuples(nEl);
    tuple[0] = 1;
    _vtkDataArray_MatId->SetTuple(0, tuple);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_MatId);
// FormulationId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_FormId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_FormId->SetName("FormulationId");
    _vtkDataArray_FormId->SetNumberOfComponents(1);
    _vtkDataArray_FormId->SetNumberOfTuples(nEl);
    tuple[0] = 900;
    _vtkDataArray_FormId->SetTuple(0, tuple);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_FormId);
// PieceId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_PieceId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_PieceId->SetName("PieceId");
    _vtkDataArray_PieceId->SetNumberOfComponents(1);
    _vtkDataArray_PieceId->SetNumberOfTuples(nEl);
    tuple[0] = 1;
    _vtkDataArray_PieceId->SetTuple(0, tuple);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PieceId);
// RegionId --------------------------------------------------------------------------------
    vtkSmartPointer<vtkIntArray> _vtkDataArray_RegionId = vtkSmartPointer<vtkIntArray>::New();
    _vtkDataArray_RegionId->SetName("RegionId");
    _vtkDataArray_RegionId->SetNumberOfComponents(1);
    _vtkDataArray_RegionId->SetNumberOfTuples(nEl);
    tuple[0] = 0;
    _vtkDataArray_RegionId->SetTuple(0, tuple);
    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_RegionId);

    unstructuredGrid->GetCellData()->AddArray(_vtkDataArray_PartId);

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
  cout << "number of subdomains: " << nSx * nSy * nSz + 1 << endl;


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

  return EXIT_SUCCESS;
}

