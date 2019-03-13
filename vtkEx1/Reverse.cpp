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

int main(int, char *[])
{
	/*vtkSmartPointer<vtkSphereSource> sphereSource =
		vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->Update();*/

	const char* lpszPathName2 = { "D:\\Project\\SurgicalGuide\\Doc\\TechnicalReview\\boolean\\SampleData\\sleeve.stl" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();

	vtkSmartPointer<vtkPolyData> mpSleevePolydata;
	mpSleevePolydata = reader->GetOutput();

	//mpSleevePolydata.Get()->

	/*double originBounds[6];
	mpSleevePolydata->GetBounds(originBounds);

	double *dflip = new double[16]{
		1.00000000,  0.00000000,  0.00000000,  0.00000000,
		0.00000000,  -1.00000000,  0.00000000,  0.00000000,
		0.00000000,  0.00000000,  -1.00000000, 15.00000000,
		0.00000000,  0.00000000,  0.00000000,  1.00000000
	};

	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
	trans->SetMatrix(dflip);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetInputData(mpSleevePolydata);
	transformFilter->SetTransform(trans);
	transformFilter->Update();*/


	vtkSTLWriter* stlWriter2 = vtkSTLWriter::New();
	const char* lpszPathName3 = { "D:\\flipSleeve1.stl" };
	stlWriter2->SetFileName(lpszPathName3);
	//stlWriter2->SetInputConnection(transformFilter->GetOutputPort());

	stlWriter2->SetInputData(mpSleevePolydata);

	stlWriter2->Write();


	return EXIT_SUCCESS;
}