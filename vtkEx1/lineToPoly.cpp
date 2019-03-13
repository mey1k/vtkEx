#include "stdafx.h"
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

int main(int, char *[])
{
	double origin[3] = { 0.0, 0.0, 0.0 };
	double p0[3] = { 1.0, 0.0, 0.0 };
	double p1[3] = { 0.0, 1.0, 0.0 };
	double p2[3] = { 0.0, 1.0, 2.0 };
	double p3[3] = { 1.0, 2.0, 3.0 };

	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(origin);
	points->InsertNextPoint(p0);
	points->InsertNextPoint(p1);
	points->InsertNextPoint(p2);
	points->InsertNextPoint(p3);

	// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> lines =
		vtkSmartPointer<vtkCellArray>::New();

	//Create four lines
	for (unsigned int i = 0; i < 4; i++)
	{
		vtkSmartPointer<vtkLine> line =
			vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0, i);
		line->GetPointIds()->SetId(1, i + 1);
		lines->InsertNextCell(line);
	}

	
	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> linesPolyData =
		vtkSmartPointer<vtkPolyData>::New();

	// Add the points to the dataset
	linesPolyData->SetPoints(points);

	// Add the lines to the dataset
	linesPolyData->SetLines(lines);

	vtkSmartPointer<vtkPolyDataMapper> inputMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	inputMapper->SetInputData(linesPolyData);

	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(1800, 900);

	// Setup both renderers
	vtkSmartPointer<vtkRenderer> leftRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(leftRenderer);
	leftRenderer->SetBackground(.6, .5, .4);

	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	leftRenderer->AddActor(inputActor);

	renderWindow->Render();
	interactor->Start();

	//return EXIT_SUCCESS;
}