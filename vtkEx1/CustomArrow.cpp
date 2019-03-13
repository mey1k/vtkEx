#define USER_MATRIX
#include "stdafx.h"
#include <vtkSmartPointer.h>

#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkArrowSource.h>
#include <vtkCellData.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkUnsignedCharArray.h>

#include <time.h>
 
static vtkSmartPointer<vtkPolyData> CreateCustomArrow();
static void ColorCellData(vtkPolyData *, vtkStdString);

int main(int, char *[])
{
  vtkSmartPointer<vtkPolyData> customArrow =
    vtkSmartPointer<vtkPolyData>::New();
  customArrow = CreateCustomArrow();

  // Generate a random start and end point
  double startPoint[3], endPoint[3];
#ifndef main
  vtkMath::RandomSeed(time(NULL));
#else
  vtkMath::RandomSeed(8775070);
#endif
  startPoint[0] = vtkMath::Random(-10,10);
  startPoint[1] = vtkMath::Random(-10,10);
  startPoint[2] = vtkMath::Random(-10,10);
  endPoint[0] = vtkMath::Random(-10,10);
  endPoint[1] = vtkMath::Random(-10,10);
  endPoint[2] = vtkMath::Random(-10,10);

  // Compute a basis
  double normalizedX[3];
  double normalizedY[3];
  double normalizedZ[3];

  // The X axis is a vector from start to end
  vtkMath::Subtract(endPoint, startPoint, normalizedX);
  double length = vtkMath::Norm(normalizedX);
  vtkMath::Normalize(normalizedX);

  // The Z axis is an arbitrary vecotr cross X
  double arbitrary[3];
  arbitrary[0] = vtkMath::Random(-10,10);
  arbitrary[1] = vtkMath::Random(-10,10);
  arbitrary[2] = vtkMath::Random(-10,10);
  vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
  vtkMath::Normalize(normalizedZ);

  // The Y axis is Z cross X
  vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();

  // Create the direction cosine matrix
  matrix->Identity();
  for (unsigned int i = 0; i < 3; i++)
    {
    matrix->SetElement(i, 0, normalizedX[i]);
    matrix->SetElement(i, 1, normalizedY[i]);
    matrix->SetElement(i, 2, normalizedZ[i]);
    }    

  // Apply the transforms
  vtkSmartPointer<vtkTransform> transform = 
    vtkSmartPointer<vtkTransform>::New();
  transform->Translate(startPoint);
  transform->Concatenate(matrix);
  transform->Scale(length, length, length);

  // Transform the polydata
  vtkSmartPointer<vtkTransformPolyDataFilter> transformPD = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformPD->SetTransform(transform);
  transformPD->SetInputData(customArrow);

  //Create a mapper and actor for the arrow
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
#ifdef USER_MATRIX
  mapper->SetInputData(customArrow);
  actor->SetUserMatrix(transform->GetMatrix());
#else
  mapper->SetInputConnection(transformPD->GetOutputPort());
#endif
  actor->SetMapper(mapper);
 
  // Create spheres for start and end point
  vtkSmartPointer<vtkSphereSource> sphereStartSource =
    vtkSmartPointer<vtkSphereSource>::New();
    sphereStartSource->SetCenter(startPoint);
  vtkSmartPointer<vtkPolyDataMapper> sphereStartMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereStartMapper->SetInputConnection(sphereStartSource->GetOutputPort());
  vtkSmartPointer<vtkActor> sphereStart =
    vtkSmartPointer<vtkActor>::New();
  sphereStart->SetMapper(sphereStartMapper);
  sphereStart->GetProperty()->SetColor(1.0, 1.0, .3);

  vtkSmartPointer<vtkSphereSource> sphereEndSource =
    vtkSmartPointer<vtkSphereSource>::New();
    sphereEndSource->SetCenter(endPoint);
  vtkSmartPointer<vtkPolyDataMapper> sphereEndMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereEndMapper->SetInputConnection(sphereEndSource->GetOutputPort());
  vtkSmartPointer<vtkActor> sphereEnd =
    vtkSmartPointer<vtkActor>::New();
  sphereEnd->SetMapper(sphereEndMapper);
  sphereEnd->GetProperty()->SetColor(1.0, .3, .3);

  //Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  //Add the actor to the scene
  renderer->AddActor(actor);
  renderer->AddActor(sphereStart);
  renderer->AddActor(sphereEnd);
  renderer->SetBackground(.1, .2, .3); // Background color dark blue

 
  //Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}

// Create a custom arrow consisting of multiple parts. All part
// positioning is done in the arrow source's coordinate system.
vtkSmartPointer<vtkPolyData> CreateCustomArrow()
{
  vtkPolyData *customArrowPd;
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  //Create an arrow part.
  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetShaftResolution(50);
  arrowSource->SetTipResolution(50); 
  arrowSource->Update();

  // Color the arrow part
  ColorCellData(arrowSource->GetOutput(), "burlywood");

  // Create sphere part
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetRadius(.06);
  sphereSource->SetPhiResolution(20);
  sphereSource->SetThetaResolution(40);
  sphereSource->Update();

  // Color the sphere part
  ColorCellData(sphereSource->GetOutput(), "burnt_sienna");

  // Position the sphere part
  vtkSmartPointer<vtkTransform> sphereTransform = 
    vtkSmartPointer<vtkTransform>::New();
  double midPoint = 1.0 / 3.0;
  sphereTransform->Translate(midPoint, 0.0, 0.0);

  // Transform the polydata
  vtkSmartPointer<vtkTransformPolyDataFilter> sphereTransformPD = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  sphereTransformPD->SetTransform(sphereTransform);
  sphereTransformPD->SetInputConnection(sphereSource->GetOutputPort());
  
  // Create the lines
  double offset = 1.0 / 3.0;
  vtkSmartPointer<vtkLineSource> line1 =
    vtkSmartPointer<vtkLineSource>::New();
  line1->SetPoint1(midPoint, - offset, 0.0);
  line1->SetPoint2(midPoint, offset, 0.0);
  line1->Update();
  ColorCellData(line1->GetOutput(), "lime");

  vtkSmartPointer<vtkLineSource> line2 =
    vtkSmartPointer<vtkLineSource>::New();
  line2->SetPoint1(midPoint, 0.0, - offset);
  line2->SetPoint2(midPoint, 0.0, offset);
  line2->Update();
  ColorCellData(line2->GetOutput(), "cyan");

  // Compose the custom arrow
  vtkSmartPointer<vtkAppendPolyData> customArrow =
    vtkSmartPointer<vtkAppendPolyData>::New();
  customArrow->AddInputConnection(arrowSource->GetOutputPort());
  customArrow->AddInputConnection(sphereTransformPD->GetOutputPort());
  customArrow->AddInputConnection(line1->GetOutputPort());
  customArrow->AddInputConnection(line2->GetOutputPort());
  customArrow->Update();

  customArrowPd = customArrow->GetOutput();

  return customArrowPd;
}

void ColorCellData(vtkPolyData *polyData, vtkStdString colorName)
{
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  // Create colors
  vtkSmartPointer<vtkUnsignedCharArray> cellData =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  cellData->SetNumberOfComponents(3);
  cellData->SetNumberOfTuples(polyData->GetNumberOfCells());
  for (int i = 0; i < polyData->GetNumberOfCells(); ++i)
    {
    unsigned char r, g, b, a;
    colors->GetColor(colorName, r, g, b, a);
    float rgb[3];
    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;
    cellData->InsertTuple(i, rgb);
    }
  polyData->GetCellData()->SetScalars(cellData);
  
}
