#include "stdafx.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkPointData.h>
#include <vtkReverseSense.h>
#include <vtkFloatArray.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkReflectionFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>

int main(int, char *[])
{
	const char* lpszPathName2 = { "D:\\Project\\SurgicalGuide\\Doc\\TechnicalReview\\boolean\\SampleData\\sleeve.stl" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();

	vtkSmartPointer<vtkPolyData> mpSleevePolydata;
	mpSleevePolydata = reader->GetOutput();
	vtkSmartPointer<vtkPoints> points = mpSleevePolydata->GetPoints();
	vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
	vtkSmartPointer<vtkDataArray> cellArray;
	vtkIdType numberofVerts = mpSleevePolydata->GetNumberOfPoints(); // vertex
	vtkIdType numberOfFaces = mpSleevePolydata->GetNumberOfCells(); // triangleCnt

	vtkIdType vertexID;

	cout << numberofVerts << endl;
	cout << numberOfFaces << endl;

	vtkSmartPointer<vtkIdList> faceIndex = vtkSmartPointer<vtkIdList>::New();

	for (int i = 0; i < numberofVerts; i++) // vertex points
	{
		double x = dataArray->GetComponent(i, 0);
		double y = dataArray->GetComponent(i, 1);
		double z = dataArray->GetComponent(i, 2);

		cout << x << ' '<< y << ' ' << z << endl;
	}

	for (int i = 0; i < numberOfFaces; i++) // cell Point ids
	{
		vtkCell* cell = mpSleevePolydata->GetCell(i);
		vtkIdList *list =  cell->GetPointIds();

		int x = list->GetId(0);
		int y = list->GetId(1);
		int z = list->GetId(2);

		cout << x << ' ' << y << ' ' << z << endl;
	}

	return EXIT_SUCCESS;
}