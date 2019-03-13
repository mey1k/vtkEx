#include "stdafx.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

#include <vtkSTLReader.h>
#include <vtkDoubleArray.h>

int main(int argc, char *argv[])
{
	vtkSmartPointer<vtkImageData> volume =
		vtkSmartPointer<vtkImageData>::New();
	double isoValue;
	/*if (argc < 3)
	{
		vtkSmartPointer<vtkSphereSource> sphereSource =
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetPhiResolution(20);
		sphereSource->SetThetaResolution(20);
		sphereSource->Update();

		double bounds[6];
		sphereSource->GetOutput()->GetBounds(bounds);
		for (unsigned int i = 0; i < 6; i += 2)
		{
			double range = bounds[i + 1] - bounds[i];
			bounds[i] = bounds[i] - .1 * range;
			bounds[i + 1] = bounds[i + 1] + .1 * range;
		}
		vtkSmartPointer<vtkVoxelModeller> voxelModeller =
			vtkSmartPointer<vtkVoxelModeller>::New();
		voxelModeller->SetSampleDimensions(50, 50, 50);
		voxelModeller->SetModelBounds(bounds);
		voxelModeller->SetScalarTypeToFloat();
		voxelModeller->SetMaximumDistance(.1);

		voxelModeller->SetInputConnection(sphereSource->GetOutputPort());
		voxelModeller->Update();
		isoValue = 0.2;
		volume->DeepCopy(voxelModeller->GetOutput());
	}
	else
	{*/

	const char* lpszPathName2 = { "D:\\target1.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();
	volume->DeepCopy(reader->GetOutput());

	vtkSmartPointer<vtkDoubleArray> weights =
		vtkSmartPointer<vtkDoubleArray>::New();
	weights->SetNumberOfValues(2);
	weights->SetValue(0, 1);
	weights->SetValue(1, 2);

	reader->GetOutput()->GetPointData()-
	reader->GetPointData()->SetScalars(weights);

	volume->DeepCopy(reader->GetOutput());
		/*vtkSmartPointer<vtkDICOMImageReader> reader =
			vtkSmartPointer<vtkDICOMImageReader>::New();
		reader->SetDirectoryName("D:\\_Data\\경북대치과 신효분 20171103\\Dicom Data");
		reader->Update();
		volume->DeepCopy(reader->GetOutput());*/
		isoValue = 0.2;
	//}

	vtkSmartPointer<vtkMarchingCubes> surface =
		vtkSmartPointer<vtkMarchingCubes>::New();

#if VTK_MAJOR_VERSION <= 5
	surface->SetInput(volume);
#else
	surface->SetInputData(volume);
#endif
	surface->ComputeNormalsOn();
	surface->SetValue(0, isoValue);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(.1, .2, .3);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(surface->GetOutputPort());
	mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	renderer->AddActor(actor);

	renderWindow->Render();
	interactor->Start();
	return EXIT_SUCCESS;
}