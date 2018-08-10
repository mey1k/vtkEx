#include "stdafx.h"

#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFillHolesFilter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkBox.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkKdTree.h>
#include <vtkKdTreePointLocator.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkImplicitModeller.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkReverseSense.h>
#include <vtkAppendPolyData.h>
#include <vtkRuledSurfaceFilter.h>
#include <vtkContourFilter.h>
#include <vtkSphereSource.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkWarpVector.h>
#include <vtkDataSetAttributes.h>
#include <vtkCleanPolyData.h>
#include <vtkStripper.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkTriangleFilter.h>
#include <vtkCamera.h>
#include <vtkDecimatePro.h>
#include <vtkParametricSpline.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkParametricFunctionSource.h>
#include <vtkFeatureEdges.h>
#include <vtkSTLWriter.h>
#include <vtkWarpScalar.h>
#include <vtkSelectPolyData.h>
#include <vtkWarpTo.h>
#include <vtkImplicitModeller.h>
#include <vtkDataArray.h>
#include <vtkVector.h>
#include <vtkVectorOperators.h>

//#include <vtkCollisionDetectionFilter>

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	bool value = 1;
	static KeyPressInteractorStyle* New();

	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
	virtual void OnKeyPress()
	{
		// Get the keypress
		std::string key = this->Interactor->GetKeySym();
		Renderer->RemoveActor(Actor);
		std::cerr << " - ";

		// "s" for "s"elect
		if (key.compare("0") == 0)
		{
			updateData->DeepCopy(inputData);
			std::cerr << "Do Initialize" << std::endl;
		}
		/*Clean Polydata*/
		else if (key.compare("1") == 0)
		{
			vtkCleanPolyData* clean =
				vtkCleanPolyData::New();
			std::cerr << "Do CleanPolyData (" << updateData->GetNumberOfPoints() << " >> ";
			clean->SetInputData(updateData);
			clean->Update();
			std::cerr << clean->GetOutput()->GetNumberOfPoints() << ")" << std::endl;


			updateData->DeepCopy(clean->GetOutput());
			clean->Delete();
		}

		/*Polydata Normals*/
		else if (key.compare("2") == 0)
		{
			// Generate normals
			vtkPolyDataNormals* normals =
				vtkPolyDataNormals::New();
			//normals->SetInputConnection(clean->GetOutputPort());
			normals->SetInputData(updateData);
			normals->SplittingOff();
			//normals->SplittingOn();
			//normals->AutoOrientNormalsOff();
			normals->ComputePointNormalsOn();
			normals->ComputeCellNormalsOn();
			//solve random normals 
			//(avoid loop traversal when VTK computes normals)
			//normals->NonManifoldTraversalOff();
			normals->ConsistencyOn();
			normals->Update();
			std::cerr << "Do vtkPolyDataNormals" << std::endl;

			updateData->DeepCopy(normals->GetOutput());
			normals->Delete();
		}

		/*ReverseSense*/
		else if (key.compare("3") == 0)
		{
			//vtkStripper* stripper =
			//	vtkStripper::New();
			//stripper->SetInputData(updateData);
			//stripper->Update();
			//std::cerr << "Do vtkStripper" << std::endl;

			//updateData->DeepCopy(stripper->GetOutput());
			//stripper->Delete();
			/*면 반전*/
			vtkReverseSense* reverseIN = vtkReverseSense::New();;
			reverseIN->SetInputData(updateData);
			//reverseIN->ReverseNormalsOn();
			reverseIN->Update();
			updateData->DeepCopy(reverseIN->GetOutput());
		}

		/*Smoothing*/
		else if (key.compare("4") == 0)
		{
			vtkSmoothPolyDataFilter* smoothFilter =
				vtkSmoothPolyDataFilter::New();
			//smoothFilter->SetInputConnection(decimate->GetOutputPort());
			smoothFilter->SetInputData(updateData);
			smoothFilter->SetNumberOfIterations(10);
			smoothFilter->SetRelaxationFactor(0.5);
			smoothFilter->FeatureEdgeSmoothingOn();
			smoothFilter->BoundarySmoothingOn();
			smoothFilter->SetFeatureAngle(30);
			//smoothFilter->GenerateErrorVectorsOn();
			smoothFilter->Update();

			std::cerr << "Do vtkSmoothPolyDataFilter( EdgeOn )  " << smoothFilter->GetFeatureAngle() << std::endl;
			updateData->DeepCopy(smoothFilter->GetOutput());
			smoothFilter->Delete();
		}
		else if (key.compare("5") == 0)
		{
			vtkFillHolesFilter* FillHoleFilter = vtkFillHolesFilter::New();
			FillHoleFilter->SetInputData(updateData);

			//FillHoleFilter->SetHoleSize(1000.0);
			FillHoleFilter->Update();

			std::cerr << "Do vtkSmoothPolyDataFilter( EdgeOff )" << std::endl;
			updateData->DeepCopy(FillHoleFilter->GetOutput());
			FillHoleFilter->Delete();

		}

		/*warp Polydata*/
		else if (key.compare("6") == 0)
		{
			vtkWarpVector* warpClipPoly =
				vtkWarpVector::New();
			warpClipPoly->SetInputData(updateData);
			warpClipPoly->SetInputArrayToProcess(0, 0, 0,
				vtkDataObject::FIELD_ASSOCIATION_POINTS,
				//vtkDataSetAttributes::VECTORS);
				vtkDataSetAttributes::NORMALS);
			warpClipPoly->SetScaleFactor(0.5);

			//warpClipPoly->SetScaleFactor(.2);
			warpClipPoly->Update();


			std::cerr << "Do vtkWarpVector( Factor .5 )" << std::endl;
			updateData->DeepCopy(warpClipPoly->GetOutput());
			warpClipPoly->Delete();
		}

		/*vtkBooleanOperationPolyDataFilter*/
		else if (key.compare("7") == 0)
		{
			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			booleanOperation->SetInputData(0, updateData);
			booleanOperation->SetInputData(1, inputData);
			booleanOperation->Update();

			std::cerr << "Do vtkBoolean( Union )" << std::endl;
			updateData->DeepCopy(booleanOperation->GetOutput());
			booleanOperation->Delete();
		}

		/*vtkDecimatePro*/
		else if (key.compare("8") == 0)
		{
			//vtkMarchingCubes* mc = vtkMarchingCubes::New();
			////mc->SetInputConnection(updateData->GetOutputPort());
			//mc->SetInputData(updateData);
			//mc->ComputeNormalsOn();
			//mc->ComputeGradientsOff();
			//mc->ComputeScalarsOff();
			//mc->SetValue(0, 3000);
			//mc->Update();
			//vtkCleanPolyData* cP = vtkCleanPolyData::New();
			//cP->SetInputData(mc->GetOutput());
			//cP->Update();

			//std::cerr << "Do vtkMarchingCubes" << std::endl;
			//updateData->DeepCopy(cP->GetOutput());
			//mc->Delete();
			//cP->Delete();
			vtkDecimatePro* decimate =
				vtkDecimatePro::New();
			decimate->SetInputData(updateData);
			//decimate->SetTargetReduction(.99); //99% reduction (if there was 100 triangles, now there will be 1)
			decimate->SetTargetReduction(.10); //10% reduction (if there was 100 triangles, now there will be 90)
			decimate->Update();
			std::cerr << "Do vtkDecimatePro" << std::endl;
			updateData->DeepCopy(decimate->GetOutput());
			decimate->Delete();
		}

		/*vtkFeatureEdges*/
		else if (key.compare("z") == 0)
		{
			vtkFeatureEdges* featureEdges = vtkFeatureEdges::New();
			//featureEdges->SetInputConnection(updateData->GetOutputPort());
			featureEdges->SetInputData(updateData);
			featureEdges->BoundaryEdgesOn();
			featureEdges->FeatureEdgesOff();
			featureEdges->ManifoldEdgesOff();
			featureEdges->NonManifoldEdgesOff();
			featureEdges->Update();
			std::cerr << "Do vtkFeatureEdges" << std::endl;
			vtkSmartPointer<vtkPolyDataMapper> outMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			outMapper->SetInputData(featureEdges->GetOutput());
			vtkSmartPointer<vtkActor> outActor =
				vtkSmartPointer<vtkActor>::New();
			outActor->SetMapper(outMapper);
			Renderer->AddActor(outActor);
			//updateData->DeepCopy(featureEdges->GetOutput());
			featureEdges->Delete();
		}
		else if (key.compare("m") == 0)
		{
			SleeveData = vtkPolyData::New();
			const char* lpszeSleevePathName = { "D:\\Sleeve.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* SleeveReader = vtkSTLReader::New();
			SleeveReader->SetFileName(lpszeSleevePathName);
			SleeveReader->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			booleanOperation->SetInputData(0, SleeveReader->GetOutput());
			booleanOperation->SetInputData(1, updateData);
			booleanOperation->Update();

			/*vtkSmartPointer<vtkPolyDataMapper> outMapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
			outMapper->SetInputConnection(booleanOperation->GetOutputPort());
			outMapper->ScalarVisibilityOff();
			vtkSmartPointer<vtkActor> outActor =
			vtkSmartPointer<vtkActor>::New();
			outActor->SetMapper(outMapper);
			Renderer->AddActor(outActor);*/
			SleeveData->DeepCopy(booleanOperation->GetOutput());
			{
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(booleanOperation->GetOutput());
				mapper->ScalarVisibilityOff();

				vtkSmartPointer<vtkProperty> backFaces =
					vtkSmartPointer<vtkProperty>::New();
				backFaces->SetSpecular(0.5);
				backFaces->SetDiffuse(.2);
				backFaces->SetAmbient(1.0);
				backFaces->SetAmbientColor(0, 0, 1);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(1., 1., 0.);
				//actor->SetBackfaceProperty(backFaces);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

				vtkSmartPointer<vtkRenderWindow> renderWindow =
					vtkSmartPointer<vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);

				vtkSmartPointer<vtkRenderWindowInteractor> interactor =
					vtkSmartPointer<vtkRenderWindowInteractor>::New();
				interactor->SetRenderWindow(renderWindow);
				renderer->AddActor(actor);

				renderer->SetBackground(.1, .2, .4);

				renderWindow->SetSize(500, 250);

				renderWindow->Render();
				interactor->Start();
			}

		}
		else if (key.compare("t") == 0)
		{
			int step = 0;
			SleeveData = vtkPolyData::New();
			const char* lpszeSleevePathName = { "D:\\Sleeve.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* SleeveReader = vtkSTLReader::New();
			SleeveReader->SetFileName(lpszeSleevePathName);
			SleeveReader->Update();

			const char* lpszeHolePathName = { "D:\\hole2.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* holereader =
				vtkSTLReader::New();
			holereader->SetFileName(lpszeHolePathName);
			holereader->Update();

			const char* lpszeOriginPathName = { "D:\\lowerjaw.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* bodyreader =
				vtkSTLReader::New();
			bodyreader->SetFileName(lpszeOriginPathName);
			bodyreader->Update();

			vtkImplicitModeller* implicitModeller =
				vtkImplicitModeller::New();
			//implicitModeller->SetSampleDimensions(50, 50, 50);
			implicitModeller->SetInputData(updateData);
			implicitModeller->AdjustBoundsOn();
			implicitModeller->SetAdjustDistance(.1); // Adjust by 10%
			implicitModeller->SetMaximumDistance(.1);
			implicitModeller->Update();

			// Compute the range to select a reasonable contour value
			double bounds[6];
			updateData->GetBounds(bounds);
			double xrange = bounds[1] - bounds[0];

			// Create the 0 isosurface
			vtkContourFilter* contourFilter = vtkContourFilter::New();
			contourFilter->SetInputConnection(implicitModeller->GetOutputPort());
			contourFilter->SetValue(0, xrange / 30.0); // 30% of xrange
			contourFilter->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);//Intersection
																							  //booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
																							  //booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, contourFilter->GetOutput());
			booleanOperation->SetInputData(1, SleeveReader->GetOutput());
			booleanOperation->Update();
			vtkReverseSense* reverseIN = vtkReverseSense::New();;
			reverseIN->SetInputData(booleanOperation->GetOutput());
			reverseIN->ReverseNormalsOn();
			reverseIN->Update();
			vtkTriangleFilter* TriangleFilter1 = vtkTriangleFilter::New();
			TriangleFilter1->SetInputData(reverseIN->GetOutput());
			TriangleFilter1->Update();
			vtkCleanPolyData* CleanPolydataA = vtkCleanPolyData::New();
			CleanPolydataA->SetInputData(TriangleFilter1->GetOutput());
			CleanPolydataA->Update();
			vtkPolyDataNormals* polydataA = vtkPolyDataNormals::New();
			polydataA->SetInputData(CleanPolydataA->GetOutput());
			polydataA->Update();
			std::cerr << "Step" << ++step << " Pass" << std::endl;

			{
				vtkSTLWriter* stlWriter = vtkSTLWriter::New();
				const char* lpszPathName = { "D:\\Result1.stl" };
				stlWriter->SetFileName(lpszPathName);
				stlWriter->SetInputConnection(polydataA->GetOutputPort());
				stlWriter->Write();
			}

			vtkBooleanOperationPolyDataFilter* booleanOperation2 =
				vtkBooleanOperationPolyDataFilter::New();
			//booleanOperation2->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//booleanOperation2->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			booleanOperation2->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation2->SetInputData(0, bodyreader->GetOutput());
			booleanOperation2->SetInputData(1, holereader->GetOutput());
			booleanOperation2->Update();

			vtkTriangleFilter* TriangleFilter2 = vtkTriangleFilter::New();
			TriangleFilter2->SetInputData(booleanOperation2->GetOutput());
			TriangleFilter2->Update();
			vtkCleanPolyData* CleanPolydataB = vtkCleanPolyData::New();
			CleanPolydataB->SetInputData(TriangleFilter2->GetOutput());
			CleanPolydataB->Update();
			vtkPolyDataNormals* polydataB = vtkPolyDataNormals::New();
			polydataB->SetInputData(CleanPolydataB->GetOutput());
			polydataB->Update();

			{
				vtkSTLWriter* stlWriter = vtkSTLWriter::New();
				const char* lpszPathName = { "D:\\Result2.stl" };
				stlWriter->SetFileName(lpszPathName);
				stlWriter->SetInputConnection(polydataB->GetOutputPort());
				stlWriter->Write();
			}
			std::cerr << "Step" << ++step << " Pass" << std::endl;
			vtkBooleanOperationPolyDataFilter* booleanOperation3 =
				vtkBooleanOperationPolyDataFilter::New();
			booleanOperation3->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//booleanOperation3->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//booleanOperation3->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			/*booleanOperation3->SetInputData(0, booleanOperation->GetOutput());
			booleanOperation3->SetInputData(1, booleanOperation2->GetOutput());*/
			booleanOperation3->SetInputConnection(0, polydataA->GetOutputPort());
			booleanOperation3->SetInputConnection(1, polydataB->GetOutputPort());
			booleanOperation3->Update();
			{
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(booleanOperation3->GetOutput());
				mapper->ScalarVisibilityOff();

				vtkSmartPointer<vtkProperty> backFaces =
					vtkSmartPointer<vtkProperty>::New();
				backFaces->SetSpecular(0.5);
				backFaces->SetDiffuse(.2);
				backFaces->SetAmbient(1.0);
				backFaces->SetAmbientColor(0, 0, 1);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(1., 1., 0.);
				//actor->SetBackfaceProperty(backFaces);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

				vtkSmartPointer<vtkRenderWindow> renderWindow =
					vtkSmartPointer<vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);

				vtkSmartPointer<vtkRenderWindowInteractor> interactor =
					vtkSmartPointer<vtkRenderWindowInteractor>::New();
				interactor->SetRenderWindow(renderWindow);
				renderer->AddActor(actor);

				renderer->SetBackground(.1, .2, .4);

				renderWindow->SetSize(500, 250);

				renderWindow->Render();
				interactor->Start();
			}
		}
		else if (key.compare("h") == 0)
		{
			const char* lpszeHolePathName = { "D:\\Result1.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* holereader =
				vtkSTLReader::New();
			holereader->SetFileName(lpszeHolePathName);
			holereader->Update();

			const char* lpszeOriginPathName = { "D:\\Result2.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* bodyreader =
				vtkSTLReader::New();
			bodyreader->SetFileName(lpszeOriginPathName);
			bodyreader->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);//Intersection
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, holereader->GetOutput());
			booleanOperation->SetInputData(1, bodyreader->GetOutput());
			booleanOperation->Update();

			{
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(booleanOperation->GetOutput());
				mapper->ScalarVisibilityOff();

				vtkSmartPointer<vtkProperty> backFaces =
					vtkSmartPointer<vtkProperty>::New();
				backFaces->SetSpecular(0.5);
				backFaces->SetDiffuse(.2);
				backFaces->SetAmbient(1.0);
				backFaces->SetAmbientColor(0, 0, 1);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(1., 1., 0.);
				//actor->SetBackfaceProperty(backFaces);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

				vtkSmartPointer<vtkRenderWindow> renderWindow =
					vtkSmartPointer<vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);

				vtkSmartPointer<vtkRenderWindowInteractor> interactor =
					vtkSmartPointer<vtkRenderWindowInteractor>::New();
				interactor->SetRenderWindow(renderWindow);
				renderer->AddActor(actor);

				renderer->SetBackground(.1, .2, .4);

				renderWindow->SetSize(500, 250);

				renderWindow->Render();
				interactor->Start();
			}
		}
		else if (key.compare("y") == 0)
		{
			const char* lpszeHolePathName = { "D:\\hole2.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* holereader =
				vtkSTLReader::New();
			holereader->SetFileName(lpszeHolePathName);
			holereader->Update();

			const char* lpszeOriginPathName = { "D:\\lowerjaw.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* bodyreader =
				vtkSTLReader::New();
			bodyreader->SetFileName(lpszeOriginPathName);
			bodyreader->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, bodyreader->GetOutput());
			booleanOperation->SetInputData(1, holereader->GetOutput());
			booleanOperation->Update();

			{
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(booleanOperation->GetOutput());
				mapper->ScalarVisibilityOff();

				vtkSmartPointer<vtkProperty> backFaces =
					vtkSmartPointer<vtkProperty>::New();
				backFaces->SetSpecular(0.5);
				backFaces->SetDiffuse(.2);
				backFaces->SetAmbient(1.0);
				backFaces->SetAmbientColor(0, 0, 1);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(1., 1., 0.);
				//actor->SetBackfaceProperty(backFaces);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

				vtkSmartPointer<vtkRenderWindow> renderWindow =
					vtkSmartPointer<vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);

				vtkSmartPointer<vtkRenderWindowInteractor> interactor =
					vtkSmartPointer<vtkRenderWindowInteractor>::New();
				interactor->SetRenderWindow(renderWindow);
				renderer->AddActor(actor);

				renderer->SetBackground(.1, .2, .4);

				renderWindow->SetSize(500, 250);

				renderWindow->Render();
				interactor->Start();
			}


		}
		else if (key.compare("a") == 0)
		{
			int step = 0;
			SleeveData = vtkPolyData::New();
			const char* lpszeSleevePathName = { "D:\\Sleeve.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* SleeveReader = vtkSTLReader::New();
			SleeveReader->SetFileName(lpszeSleevePathName);
			SleeveReader->Update();

			std::cerr << "@ " << step++ << std::endl;

			vtkBooleanOperationPolyDataFilter* booleanOperation1 =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation1->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			booleanOperation1->SetInputData(0, SleeveReader->GetOutput());
			booleanOperation1->SetInputData(1, updateData);
			booleanOperation1->Update();
			SleeveData->DeepCopy(booleanOperation1->GetOutput());
			std::cerr << "@ " << step++ << std::endl;

			const char* lpszeHolePathName = { "D:\\hole2.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* holereader =
				vtkSTLReader::New();
			holereader->SetFileName(lpszeHolePathName);
			holereader->Update();
			std::cerr << "@ " << step++ << std::endl;
			vtkImplicitModeller* implicitModeller =
				vtkImplicitModeller::New();
			//implicitModeller->SetSampleDimensions(50, 50, 50);
			implicitModeller->SetInputData(updateData);
			implicitModeller->AdjustBoundsOn();
			implicitModeller->SetAdjustDistance(.1); // Adjust by 10%
			implicitModeller->SetMaximumDistance(.1);
			implicitModeller->Update();
			std::cerr << "@ " << step++ << std::endl;
			// Compute the range to select a reasonable contour value
			double bounds[6];
			updateData->GetBounds(bounds);
			double xrange = bounds[1] - bounds[0];

			// Create the 0 isosurface
			vtkContourFilter* contourFilter = vtkContourFilter::New();
			contourFilter->SetInputConnection(implicitModeller->GetOutputPort());
			contourFilter->SetValue(0, xrange / 30.0); // 30% of xrange
			contourFilter->Update();
			std::cerr << "@ " << step++ << std::endl;
			vtkPolyData* result = vtkPolyData::New();
			result->DeepCopy(contourFilter->GetOutput());
			std::cerr << "@ " << step++ << std::endl;
			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, contourFilter->GetOutput());
			booleanOperation->SetInputData(1, holereader->GetOutput());
			booleanOperation->Update();
			std::cerr << "@ " << step++ << std::endl;
			vtkReverseSense* reverseIN = vtkReverseSense::New();;
			reverseIN->SetInputData(booleanOperation->GetOutput());
			reverseIN->ReverseNormalsOn();
			reverseIN->Update();
			std::cerr << "@ " << step++ << std::endl;

			//vtkAppendPolyData* append = vtkAppendPolyData::New();
			////append->AddInputConnection(booleanOperation->GetOutputPort());
			//append->AddInputData(reverseIN->GetOutput());
			//append->AddInputData(SleeveData);
			//append->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation2 =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//booleanOperation2->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation2->SetInputData(0, reverseIN->GetOutput());
			booleanOperation2->SetInputData(1, SleeveData);
			booleanOperation2->Update();


			vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			const char* lpszPathName = { "D:\\Result1.stl" };
			stlWriter->SetFileName(lpszPathName);
			stlWriter->SetInputConnection(booleanOperation2->GetOutputPort());
			stlWriter->Write();

			std::cerr << "@ " << step++ << std::endl;
			std::cerr << "Do vtkFeatureEdges+warp+Hole Boolean" << std::endl;
			updateData->DeepCopy(booleanOperation2->GetOutput());

		}
		else if (key.compare("x") == 0)
		{
			const char* lpszeHolePathName = { "D:\\hole2.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* holereader =
				vtkSTLReader::New();
			holereader->SetFileName(lpszeHolePathName);
			holereader->Update();

			vtkImplicitModeller* implicitModeller =
				vtkImplicitModeller::New();
			//implicitModeller->SetSampleDimensions(50, 50, 50);
			implicitModeller->SetInputData(updateData);
			implicitModeller->AdjustBoundsOn();
			implicitModeller->SetAdjustDistance(.1); // Adjust by 10%
			implicitModeller->SetMaximumDistance(.1);
			implicitModeller->Update();

			// Compute the range to select a reasonable contour value
			double bounds[6];
			updateData->GetBounds(bounds);
			double xrange = bounds[1] - bounds[0];

			// Create the 0 isosurface
			vtkContourFilter* contourFilter = vtkContourFilter::New();
			contourFilter->SetInputConnection(implicitModeller->GetOutputPort());
			contourFilter->SetValue(0, xrange / 30.0); // 30% of xrange
			contourFilter->Update();

			vtkPolyData* result = vtkPolyData::New();
			result->DeepCopy(contourFilter->GetOutput());

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, contourFilter->GetOutput());
			booleanOperation->SetInputData(1, holereader->GetOutput());
			booleanOperation->Update();

			updateData->DeepCopy(booleanOperation->GetOutput());
			{
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputData(booleanOperation->GetOutput());
				mapper->ScalarVisibilityOff();

				vtkSmartPointer<vtkProperty> backFaces =
					vtkSmartPointer<vtkProperty>::New();
				backFaces->SetSpecular(0.5);
				backFaces->SetDiffuse(.2);
				backFaces->SetAmbient(1.0);
				backFaces->SetAmbientColor(0, 0, 1);

				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(1., 1., 0.);
				//actor->SetBackfaceProperty(backFaces);

				vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

				vtkSmartPointer<vtkRenderWindow> renderWindow =
					vtkSmartPointer<vtkRenderWindow>::New();
				renderWindow->AddRenderer(renderer);

				vtkSmartPointer<vtkRenderWindowInteractor> interactor =
					vtkSmartPointer<vtkRenderWindowInteractor>::New();
				interactor->SetRenderWindow(renderWindow);
				renderer->AddActor(actor);

				renderer->SetBackground(.1, .2, .4);

				renderWindow->SetSize(500, 250);

				renderWindow->Render();
				interactor->Start();
			}
			std::cerr << "Do vtkFeatureEdges+warp+Hole Boolean" << std::endl;

			//vtkSmartPointer<vtkPolyDataMapper> outMapper =
			//vtkSmartPointer<vtkPolyDataMapper>::New();
			//outMapper->SetInputConnection(booleanOperation->GetOutputPort());
			//outMapper->ScalarVisibilityOff();
			//vtkSmartPointer<vtkActor> outActor =
			//vtkSmartPointer<vtkActor>::New(); 
			//outActor->SetMapper(outMapper);
			//Renderer->AddActor(outActor);
			//updateData->DeepCopy(featureEdges->GetOutput());
			//holereader->Delete();
			/*	holereader->Delete();
			implicitModeller->Delete();
			result->Delete();
			contourFilter->Delete();
			booleanOperation->Delete();*/
			return;
		}
		else if (key.compare("c") == 0)
		{
			//vtkAppendPolyData* append = vtkAppendPolyData::New();
			////append->AddInputConnection(booleanOperation->GetOutputPort());
			//append->AddInputData(updateData);
			//append->AddInputData(SleeveData);
			//append->Update();

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, updateData);
			booleanOperation->SetInputData(1, SleeveData);
			booleanOperation->Update();

			//vtkMarchingCubes* pMarchingCube = vtkMarchingCubes::New();
			//pMarchingCube->SetInputData(append->GetOutput());
			//pMarchingCube->SetValue(0, 2020);
			//pMarchingCube->ComputeGradientsOff();
			//pMarchingCube->ComputeNormalsOn();	// Normal을 계산하면 rendering 결과 좋아짐
			//pMarchingCube->ComputeScalarsOff();
			//pMarchingCube->Update();
			//{
			//	vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//	const char* lpszPathName = { "D:\\_Append.stl" };
			//	stlWriter->SetFileName(lpszPathName);
			//	stlWriter->SetInputConnection(append->GetOutputPort());
			//	stlWriter->Write();
			//}
			//{
			//	vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//	const char* lpszPathName = { "D:\\_Body.stl" };
			//	stlWriter->SetFileName(lpszPathName);
			//	stlWriter->SetInputData(updateData);
			//	stlWriter->Write();
			//}
			//{
			//	vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//	const char* lpszPathName = { "D:\\_Sleeve.stl" };
			//	stlWriter->SetFileName(lpszPathName);
			//	stlWriter->SetInputData(SleeveData);
			//	stlWriter->Write();
			//}

			updateData->DeepCopy(booleanOperation->GetOutput());

		}
		else if (key.compare("p") == 0)
		{
			/*Transform 위치 설정 및 적용*/
			vtkTransform* xfm = vtkTransform::New();;
			xfm->Translate(0, 0, -20);/*나중에 CutPlane Vector 값으로 변경 필요*/
			vtkTransformPolyDataFilter* xfmPdMax = vtkTransformPolyDataFilter::New();;
			xfmPdMax->SetInputData(inputData);
			xfmPdMax->SetTransform(xfm);
			xfmPdMax->Update();

			/*면 반전*/
			vtkReverseSense* reverseIN = vtkReverseSense::New();;
			reverseIN->SetInputConnection(xfmPdMax->GetOutputPort());
			//reverseIN->ReverseNormalsOn();
			reverseIN->Update();

			/*Transform 적용된 데이터의 outline 확보*/
			vtkFeatureEdges* featureEdgesTrans = vtkFeatureEdges::New();
			featureEdgesTrans->SetInputData(reverseIN->GetOutput());
			featureEdgesTrans->BoundaryEdgesOn();
			featureEdgesTrans->FeatureEdgesOff();
			featureEdgesTrans->ManifoldEdgesOff();
			featureEdgesTrans->NonManifoldEdgesOff();
			featureEdgesTrans->Update();

			// Compute the range to select a reasonable contour value
			double bounds[6];
			featureEdgesTrans->GetOutput()->GetBounds(bounds);
			double xrange = bounds[1] - bounds[0];
			double yrange = bounds[3] - bounds[2];
			double zrange = bounds[5] - bounds[4];

			{
				vtkSphereSource* sphereSource =
					vtkSphereSource::New();
				sphereSource->SetCenter(bounds[1] / 2, bounds[3] / 2, bounds[5] / 2);
				sphereSource->Update();
				vtkPolyDataMapper*outMapper =
					vtkPolyDataMapper::New();
				outMapper->SetInputConnection(sphereSource->GetOutputPort());
				outMapper->ScalarVisibilityOff();
				vtkActor* outActor =
					vtkActor::New();
				outActor->SetMapper(outMapper);
				outActor->GetProperty()->SetColor(1, 0, 0);

				Renderer->AddActor(outActor);
				Renderer->GetRenderWindow()->Render();

				sphereSource->Delete();
				outMapper->Delete();
				outActor->Delete();
			}

			{
				vtkSphereSource* sphereSource =
					vtkSphereSource::New();
				sphereSource->SetCenter(bounds[0], bounds[2], bounds[4]);
				sphereSource->Update();
				vtkPolyDataMapper*outMapper =
					vtkPolyDataMapper::New();
				outMapper->SetInputConnection(sphereSource->GetOutputPort());
				outMapper->ScalarVisibilityOff();
				vtkActor* outActor =
					vtkActor::New();
				outActor->SetMapper(outMapper);
				outActor->GetProperty()->SetColor(0, 1, 0);

				Renderer->AddActor(outActor);
				Renderer->GetRenderWindow()->Render();

				sphereSource->Delete();
				outMapper->Delete();
				outActor->Delete();
			}

			vtkPolyDataConnectivityFilter* TransEdgeConFilter = vtkPolyDataConnectivityFilter::New();
			TransEdgeConFilter->SetInputConnection(featureEdgesTrans->GetOutputPort());
			TransEdgeConFilter->SetExtractionModeToLargestRegion();
			TransEdgeConFilter->Update();

			vtkFeatureEdges* featureEdges = vtkFeatureEdges::New();
			featureEdges->SetInputData(inputData);
			featureEdges->BoundaryEdgesOn();
			featureEdges->FeatureEdgesOff();
			featureEdges->ManifoldEdgesOff();
			featureEdges->NonManifoldEdgesOff();
			featureEdges->Update();

			vtkPolyDataConnectivityFilter* confilter = vtkPolyDataConnectivityFilter::New();
			confilter->SetInputConnection(featureEdges->GetOutputPort());
			confilter->SetExtractionModeToLargestRegion();
			confilter->Update();

			std::cerr << "@featureEdges : " << featureEdges->GetOutput()->GetNumberOfPoints() << std::endl;
			std::cerr << "@featureEdgesTrans : " << featureEdgesTrans->GetOutput()->GetNumberOfPoints() << std::endl;
			//vtkKdTreePointLocator* pointTreeEdge = vtkKdTreePointLocator::New();
			//vtkKdTreePointLocator* pointTreeTransEdge = vtkKdTreePointLocator::New();
			////pointTreeEdge->BuildLocatorFromPoints(featureEdges->GetOutput()->GetPoints());
			//pointTreeEdge->SetDataSet(featureEdges->GetOutput());
			//pointTreeEdge->BuildLocator();
			////pointTreeTransEdge->BuildLocatorFromPoints(featureEdgesTrans->GetOutput()->GetPoints());
			//pointTreeTransEdge->SetDataSet(featureEdgesTrans->GetOutput());
			//pointTreeTransEdge->BuildLocator();
			vtkLinearExtrusionFilter* extrude =
				vtkLinearExtrusionFilter::New();
			extrude->SetInputConnection(confilter->GetOutputPort());
			//extrude->SetExtrusionTypeToNormalExtrusion();
			extrude->SetVector(0, 0, -20);
			//extrude->SetScaleFactor(0.5);
			extrude->Update();

			vtkPolyDataNormals* Exnormals = vtkPolyDataNormals::New();
			Exnormals->SetInputData(extrude->GetOutput());
			Exnormals->SplittingOff();
			Exnormals->ComputePointNormalsOn();
			Exnormals->ComputeCellNormalsOn();
			Exnormals->ConsistencyOn();
			Exnormals->Update();

			vtkReverseSense* reverseEx = vtkReverseSense::New();;
			reverseEx->SetInputConnection(Exnormals->GetOutputPort());
			//reverseIN->ReverseNormalsOn();
			reverseEx->Update();

			////vtkIdList* result = vtkIdList::New();
			//vtkIdType result;
			//double p[3], preP[3];
			//vtkPoints* OriginPoints = vtkPoints::New();
			//for (int i = 0; i < featureEdges->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
			//{
			//	if (i == 0)
			//	{
			//		//pointTreeEdge->FindClosestNPoints(1, featureEdges->GetOutput()->GetPoints()->GetPoint(i), result);
			//		result = pointTreeEdge->FindClosestPoint(featureEdges->GetOutput()->GetPoints()->GetPoint(i));
			//		featureEdges->GetOutput()->GetPoints()->GetPoint(result, p);
			//		ReallyDeletePoint(featureEdges->GetOutput()->GetPoints(), result);
			//	}
			//	else
			//	{
			//		//pointTreeEdge->FindClosestNPoints(1, p, result);
			//		result = pointTreeEdge->FindClosestPoint(p);
			//		featureEdges->GetOutput()->GetPoints()->GetPoint(result, p);
			//		ReallyDeletePoint(featureEdges->GetOutput()->GetPoints(), result);
			//	}
			//	{
			//		vtkSphereSource* sphereSource =
			//		vtkSphereSource::New();
			//		sphereSource->SetCenter(p);
			//		sphereSource->Update();
			//		vtkPolyDataMapper*outMapper =
			//		vtkPolyDataMapper::New();
			//		outMapper->SetInputConnection(sphereSource->GetOutputPort());
			//		outMapper->ScalarVisibilityOff();
			//		vtkActor* outActor =
			//		vtkActor::New();
			//		outActor->SetMapper(outMapper);
			//		outActor->GetProperty()->SetColor(1, 0, 0);

			//		Renderer->AddActor(outActor);
			//		Renderer->GetRenderWindow()->Render();
			//		
			//		sphereSource->Delete();
			//		outMapper->Delete();
			//		outActor->Delete();
			//	}
			//	OriginPoints->InsertNextPoint(p);
			//	std::cerr << "p " << i << " : " << p[0] << "  " << p[1] << "  " << p[2] << std::endl;

			//	//vtkIdType Warp_point_id = result->GetId(0);
			//}

			//vtkPoints* TransPoints = vtkPoints::New();
			//for (int i = 0; i < featureEdgesTrans->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
			//{
			//	if (i == 0)
			//	{
			//		result = pointTreeTransEdge->FindClosestPoint(featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(i));
			//		featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(result, p);
			//		ReallyDeletePoint(featureEdgesTrans->GetOutput()->GetPoints(), result);
			//	}
			//	else
			//	{
			//		result = pointTreeTransEdge->FindClosestPoint(p);
			//		featureEdges->GetOutput()->GetPoints()->GetPoint(result, p);
			//		featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(result, p);
			//			/*featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(result->GetId(1), p);
			//			std::cerr << "p 1 : " << p[0] << "  " << p[1] << "  " << p[2] << std::endl;
			//			featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(result->GetId(2), p);
			//			std::cerr << "p 2 : " << p[0] << "  " << p[1] << "  " << p[2] << std::endl << std::endl;*/

			//	}
			//	TransPoints->InsertNextPoint(p);
			//	std::cerr << "p " << i << " : " << p[0] << "  " << p[1] << "  " << p[2] << std::endl;

			//	//std::cerr << "p : " << p[0] <<"  " << p[1] << "  " << p[2] << std::endl;
			//	preP[0] = p[0];
			//	preP[1] = p[1];
			//	preP[2] = p[2];
			//	//vtkIdType Warp_point_id = result->GetId(0);
			//}
			/*vtkTriangleFilter* TriangleF = vtkTriangleFilter::New();
			TriangleF->SetInputConnection(TransEdgeConFilter->GetOutputPort());
			TriangleF->Update();*/
			vtkAppendPolyData* appendPD = vtkAppendPolyData::New();
			appendPD->AddInputData(inputData);
			appendPD->AddInputConnection(reverseIN->GetOutputPort());
			appendPD->AddInputConnection(reverseEx->GetOutputPort());
			appendPD->Update();

			vtkFillHolesFilter* FillHoleFilter = vtkFillHolesFilter::New();
			FillHoleFilter->SetInputConnection(appendPD->GetOutputPort());
			FillHoleFilter->SetHoleSize(1000.0);
			FillHoleFilter->Update();
			vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			const char* lpszPathName = { "D:\\TestBooleansTargetResult.stl" };
			stlWriter->SetFileName(lpszPathName);
			stlWriter->SetInputConnection(FillHoleFilter->GetOutputPort());
			stlWriter->Write();

			//vtkKdTree* pointTreeWarpPolyData = vtkKdTree::New();
			//pointTreeWarpPolyData->BuildLocatorFromPoints(inputData->GetPoints());

			//vtkKdTree* pointTreeOriginPolyData = vtkKdTree::New();
			//pointTreeOriginPolyData->BuildLocatorFromPoints(reverseIN->GetOutput()->GetPoints());

			//vtkKdTree* pointTree_AppendPolyData = vtkKdTree::New();
			//pointTree_AppendPolyData->BuildLocatorFromPoints(appendPD->GetOutput()->GetPoints());

			//vtkIdList* resultList = vtkIdList::New();
			//vtkCellArray* Polys = vtkCellArray::New();
			//Polys->DeepCopy(appendPD->GetOutput()->GetPolys());
			//vtkPoints* ClosestWarpPolyDataPoints = vtkPoints::New();
			//vtkPoints* ClosestOriginPolyDataPoints = vtkPoints::New();
			////double p[3];

			////for (int i = 0; i < featureEdges->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
			////{
			////	pointTreeWarpPolyData->FindClosestNPoints(1, featureEdges->GetOutput()->GetPoints()->GetPoint(i), result);
			////	inputData->GetPoints()->GetPoint(result->GetId(0), p);
			////	ClosestWarpPolyDataPoints->InsertNextPoint(p);
			////}
			////for (int i = 0; i <featureEdgesTrans->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
			////{
			////	pointTreeOriginPolyData->FindClosestNPoints(1, featureEdgesTrans->GetOutput()->GetPoints()->GetPoint(i), result);
			////	reverseIN->GetOutput()->GetPoints()->GetPoint(result->GetId(0), p);
			////	ClosestOriginPolyDataPoints->InsertNextPoint(p);
			////}

			//vtkIdType Pre_WarpPointID;
			//vtkIdType Pre_OriginPointID;
			//vtkCellArray* LineArray = vtkCellArray::New();
			//for (int i = 0; i <featureEdges->GetOutput()->GetPoints()->GetNumberOfPoints()-1; i++)//10; i++)//
			//{
			//	pointTree_AppendPolyData->FindClosestNPoints(1, OriginPoints->GetPoint(i), resultList);
			//	vtkIdType Warp_point_id = resultList->GetId(0);
			//	pointTree_AppendPolyData->FindClosestNPoints(1, TransPoints->GetPoint(i), resultList);
			//	vtkIdType Origin_point_id = resultList->GetId(0);
			//	std::cerr << "Origin_point_id : " << OriginPoints->GetPoint(i)[0] << OriginPoints->GetPoint(i)[1] << OriginPoints->GetPoint(i)[2] << std::endl;
			//	std::cerr << "Warp_point_id : " << Warp_point_id << std::endl << std::endl;
			//	if (i > 0)
			//	{
			//	//	vtkTriangle* Triangle1 = vtkTriangle::New();
			//	//	Triangle1->GetPointIds()->SetId(0, Origin_point_id);
			//	//	Triangle1->GetPointIds()->SetId(1, Pre_OriginPointID);
			//	//	Triangle1->GetPointIds()->SetId(2, Pre_WarpPointID);
			//	//	Polys->InsertNextCell(Triangle1);
			//	//	vtkTriangle* Triangle2 = vtkTriangle::New();
			//	//	Triangle2->GetPointIds()->SetId(0, Pre_WarpPointID);
			//	//	Triangle2->GetPointIds()->SetId(1, Warp_point_id);
			//	//	Triangle2->GetPointIds()->SetId(2, Origin_point_id);
			//	//	Polys->InsertNextCell(Triangle2);
			//	}
			//	//else
			//	//{
			//		if (i > 0)
			//		{
			//			vtkLine* Line1 = vtkLine::New();
			//			Line1->GetPointIds()->SetId(0, Pre_WarpPointID);
			//			Line1->GetPointIds()->SetId(1, Origin_point_id);
			//			vtkLine* Line2 = vtkLine::New();
			//			Line1->GetPointIds()->SetId(0, Origin_point_id);
			//			Line1->GetPointIds()->SetId(1, Pre_OriginPointID);
			//			Polys->InsertNextCell(Line1);
			//			Polys->InsertNextCell(Line2);
			//		}
			//		vtkLine* Line = vtkLine::New();
			//		Line->GetPointIds()->SetId(0, Origin_point_id);
			//		Line->GetPointIds()->SetId(1, Warp_point_id);
			//		Polys->InsertNextCell(Line);
			//	//}
			//	Pre_WarpPointID = Warp_point_id;
			//	Pre_OriginPointID = Origin_point_id;
			//}
			//appendPD->GetOutput()->SetPolys(Polys);
			//appendPD->Update();

			//
			//// Merge two lines on one polydata 
			//vtkAppendPolyData * mergedPolyData = vtkAppendPolyData::New();
			//mergedPolyData->AddInputConnection(featureEdges->GetOutputPort());
			//mergedPolyData->AddInputConnection(featureEdgesTrans->GetOutputPort());
			//mergedPolyData->Update();

			/*		vtkRuledSurfaceFilter* ruledSurfaceFilter = vtkRuledSurfaceFilter::New();
			ruledSurfaceFilter->SetInputData(mergedPolyData->GetOutput());
			ruledSurfaceFilter->SetResolution(featureEdges->GetOutput()->GetPoints()->GetNumberOfPoints(),20);
			ruledSurfaceFilter->SetRuledModeToResample();
			ruledSurfaceFilter->Update();*/

			/*Start 확인 (Sphere) */
			//for (int i = 0; i < 5; i++)
			//{
			//	vtkSmartPointer<vtkSphereSource> sphereSource =
			//	vtkSmartPointer<vtkSphereSource>::New();
			//	sphereSource->SetCenter(featureEdges->GetOutput()->GetPoint(i));
			//	sphereSource->Update();
			//	vtkSmartPointer<vtkSphereSource> sphereSource2 =
			//	vtkSmartPointer<vtkSphereSource>::New();
			//	sphereSource2->SetCenter(featureEdgesTrans->GetOutput()->GetPoint(i));
			//	sphereSource2->Update();
			//	{
			//	vtkSmartPointer<vtkPolyDataMapper> outMapper =
			//	vtkSmartPointer<vtkPolyDataMapper>::New();
			//	outMapper->SetInputConnection(sphereSource2->GetOutputPort());
			//	outMapper->ScalarVisibilityOff();
			//	vtkSmartPointer<vtkActor> outActor =
			//	vtkSmartPointer<vtkActor>::New();
			//	outActor->SetMapper(outMapper);
			//	outActor->GetProperty()->SetColor(1, 0, 0);

			//	Renderer->AddActor(outActor);
			//	}
			//	{
			//	vtkSmartPointer<vtkPolyDataMapper> outMapper =
			//	vtkSmartPointer<vtkPolyDataMapper>::New();
			//	outMapper->SetInputConnection(sphereSource->GetOutputPort());
			//	outMapper->ScalarVisibilityOff();
			//	vtkSmartPointer<vtkActor> outActor =
			//	vtkSmartPointer<vtkActor>::New();
			//	outActor->SetMapper(outMapper);
			//	outActor->GetProperty()->SetColor(1, 0, 0);

			//	Renderer->AddActor(outActor);
			//	}
			//}


			//{
			//	vtkSmartPointer<vtkPolyDataMapper> outMapper =
			//		vtkSmartPointer<vtkPolyDataMapper>::New();
			//	outMapper->SetInputConnection(featureEdgesTrans->GetOutputPort());
			//	outMapper->ScalarVisibilityOff();
			//	vtkSmartPointer<vtkActor> outActor =
			//		vtkSmartPointer<vtkActor>::New();
			//	outActor->SetMapper(outMapper);
			//	outActor->GetProperty()->SetColor(1, 0, 0);

			//	Renderer->AddActor(outActor);
			//}
			//{
			//	vtkSmartPointer<vtkPolyDataMapper> outMapper =
			//		vtkSmartPointer<vtkPolyDataMapper>::New();
			//	outMapper->SetInputConnection(featureEdges->GetOutputPort());
			//	outMapper->ScalarVisibilityOff();
			//	vtkSmartPointer<vtkActor> outActor =
			//		vtkSmartPointer<vtkActor>::New();
			//	outActor->SetMapper(outMapper);
			//	outActor->GetProperty()->SetColor(0, 0, 1);

			//	Renderer->AddActor(outActor);
			//}


			//vtkPolyDataNormals* normals = 	vtkPolyDataNormals::New();
			//normals->SetInputData(inputData);
			//normals->SplittingOff();
			//normals->ComputePointNormalsOn();
			//normals->ComputeCellNormalsOn();
			//normals->ConsistencyOn();
			//normals->Update();

			//vtkWarpVector * warpClipPoly = vtkWarpVector::New();
			//warpClipPoly->SetInputData(normals->GetOutput());
			//warpClipPoly->SetInputArrayToProcess(0, 0, 0,
			//	vtkDataObject::FIELD_ASSOCIATION_POINTS,
			//	vtkDataSetAttributes::NORMALS);
			//warpClipPoly->SetScaleFactor(0.01);
			//warpClipPoly->Update();
			{
				//vtkSmartPointer<vtkPolyDataMapper> outMapper =
				//	vtkSmartPointer<vtkPolyDataMapper>::New();
				//outMapper->SetInputConnection(appendPD->GetOutputPort());
				//outMapper->ScalarVisibilityOff();
				//vtkSmartPointer<vtkActor> outActor =
				//	vtkSmartPointer<vtkActor>::New();
				//outActor->SetMapper(outMapper);
				//outActor->GetProperty()->SetColor(1, 1, 0);
				//vtkSmartPointer<vtkProperty> backFaces1 =
				//	vtkSmartPointer<vtkProperty>::New();
				//backFaces1->SetSpecular(0.5);
				//backFaces1->SetDiffuse(.2);
				//backFaces1->SetAmbient(1.0);
				//backFaces1->SetAmbientColor(.0, .0, .0);
				//backFaces1->SetColor(.3, .3, .3);
				//outActor->SetBackfaceProperty(backFaces1);
				//Renderer->AddActor(outActor);
			}
			/*{
			vtkSmartPointer<vtkPolyDataMapper> outMapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
			outMapper->SetInputConnection(reverseIN->GetOutputPort());
			outMapper->ScalarVisibilityOff();
			vtkSmartPointer<vtkActor> outActor =
			vtkSmartPointer<vtkActor>::New();
			outActor->SetMapper(outMapper);
			outActor->GetProperty()->SetColor(1, 1, 0);
			vtkSmartPointer<vtkProperty> backFaces1 =
			vtkSmartPointer<vtkProperty>::New();
			backFaces1->SetSpecular(0.5);
			backFaces1->SetDiffuse(.2);
			backFaces1->SetAmbient(1.0);
			backFaces1->SetAmbientColor(.0, .0, .0);
			backFaces1->SetColor(.3, .3, .3);
			outActor->SetBackfaceProperty(backFaces1);
			Renderer->AddActor(outActor);
			}*/
			BooleanTargetData = vtkPolyData::New();

			std::cerr << "Do Make Boolean Source" << std::endl;
			BooleanTargetData->DeepCopy(FillHoleFilter->GetOutput());

			xfm->Delete();
			xfmPdMax->Delete();
			reverseIN->Delete();
			featureEdges->Delete();
			featureEdgesTrans->Delete();
			//ruledSurfaceFilter->Delete();
			/*normals->Delete();
			warpClipPoly->Delete();*/
			appendPD->Delete();
			stlWriter->Delete();
			extrude->Delete();
		}
		else if (key.compare("f") == 0)
		{
			vtkFillHolesFilter* FillHoleFilter = vtkFillHolesFilter::New();
			FillHoleFilter->SetInputData(updateData);

			//FillHoleFilter->SetHoleSize(1000.0);
			FillHoleFilter->Update();

			updateData->DeepCopy(FillHoleFilter->GetOutput());
			FillHoleFilter->Delete();
		}
		else if (key.compare("v") == 0)
		{
			int step = 0;
			const char* lpszeOriginPathName = { "D:\\lowerjaw.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* bodyreader =
				vtkSTLReader::New();
			bodyreader->SetFileName(lpszeOriginPathName);
			bodyreader->Update();
			std::cerr << "@ " << step++ << std::endl;
			double bounds[6];
			updateData->GetBounds(bounds);
			std::cerr << "@ " << step++ << std::endl;
			//vtkBox* implicitCube =
			//	vtkBox::New();
			//implicitCube->SetBounds(bounds);

			/*	vtkTransform* trans = vtkTransform::New();
			trans->Scale(1.2, 1.2, 1.2);

			vtkTransformPolyDataFilter* TransFilter = vtkTransformPolyDataFilter::New();
			TransFilter->SetInputDataObject(implicitCube);*/

			//vtkClipPolyData* KnifeBoundsClip = vtkClipPolyData::New();
			//KnifeBoundsClip->SetInputData(holereader->GetOutput());
			//KnifeBoundsClip->SetClipFunction(implicitCube);
			//KnifeBoundsClip->InsideOutOn();
			//KnifeBoundsClip->Update();

			//{

			//	vtkBooleanOperationPolyDataFilter* booleanOperation =
			//		vtkBooleanOperationPolyDataFilter::New();
			//	//	booleanOperation->SetOperationToUnion();
			//	//	booleanOperation->SetOperationToIntersection();
			//	booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//	booleanOperation->SetInputData(0, updateData);
			//	booleanOperation->SetInputData(1, BooleanTargetData);
			//	booleanOperation->Update();
			//	vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//	const char* lpszPathName = { "D:\\TestResult(intersection).stl" };
			//	stlWriter->SetFileName(lpszPathName);
			//	stlWriter->SetInputConnection(booleanOperation->GetOutputPort());
			//	stlWriter->Write();
			//	booleanOperation->Delete();
			//}
			//{
			//	vtkBooleanOperationPolyDataFilter* booleanOperation =
			//		vtkBooleanOperationPolyDataFilter::New();
			//	//	booleanOperation->SetOperationToUnion();
			//	//	booleanOperation->SetOperationToIntersection();
			//	booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//	booleanOperation->SetInputData(0, updateData);v
			//	booleanOperation->SetInputData(1, BooleanTargetData);
			//	booleanOperation->Update();
			//	vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//	const char* lpszPathName = { "D:\\TestResult(difference).stl" };
			//	stlWriter->SetFileName(lpszPathName);
			//	stlWriter->SetInputConnection(booleanOperation->GetOutputPort());
			//	stlWriter->Write();
			//	booleanOperation->Delete();
			//}

			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(0, updateData);
			booleanOperation->SetInputData(1, bodyreader->GetOutput());
			booleanOperation->Update();

			//vtkPolyDataConnectivityFilter* confilter = vtkPolyDataConnectivityFilter::New();
			//confilter->SetInputConnection(booleanOperation->GetOutputPort());
			//confilter->SetExtractionModeToLargestRegion();
			//confilter->Update();

			//vtkFillHolesFilter* FillHoleFilter = vtkFillHolesFilter::New();
			//FillHoleFilter->SetInputData(confilter->GetOutput());
			//FillHoleFilter->Update();

			std::cerr << "@ " << step++ << std::endl;

			vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			const char* lpszPathName = { "D:\\TestResult(union).stl" };
			stlWriter->SetFileName(lpszPathName);
			stlWriter->SetInputConnection(booleanOperation->GetOutputPort());
			stlWriter->Write();

			std::cerr << "@ " << step++ << std::endl;

			updateData->DeepCopy(booleanOperation->GetOutput());
			booleanOperation->Delete();
			//confilter->Delete();

			std::cerr << "Do vtkFeatureEdges+warp+Hole Boolean" << std::endl;
			//booleanOperation->Delete();
		}
		else if (key.compare("A") == 0)
		{
			std::cerr << "Sleeve + Body Boolean" << std::endl;
			const char* lpszeSleevePathName = { "D:\\_Sleeve.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* SleeveReader =
				vtkSTLReader::New();
			SleeveReader->SetFileName(lpszeSleevePathName);
			SleeveReader->Update();

			const char* lpszeBodyPathName = { "D:\\_Body.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
			vtkSTLReader* BodyReader =
				vtkSTLReader::New();
			BodyReader->SetFileName(lpszeBodyPathName);
			BodyReader->Update();

			//vtkAppendPolyData* append = vtkAppendPolyData::New();
			////append->AddInputConnection(booleanOperation->GetOutputPort());
			//append->AddInputData(SleeveReader->GetOutput());
			//append->AddInputData(BodyReader->GetOutput());
			//append->Update();

			//vtkCleanPolyData* cleanFilter =
			//	vtkCleanPolyData::New();
			//cleanFilter->SetInputConnection(append->GetOutputPort());
			//cleanFilter->Update();


			vtkBooleanOperationPolyDataFilter* booleanOperation =
				vtkBooleanOperationPolyDataFilter::New();
			//	booleanOperation->SetOperationToUnion();
			//	booleanOperation->SetOperationToIntersection();
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
			booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION);
			//booleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
			booleanOperation->SetInputData(1, BodyReader->GetOutput());
			booleanOperation->SetInputData(0, SleeveReader->GetOutput());
			booleanOperation->Update();

			updateData->DeepCopy(booleanOperation->GetOutput());

		}
		else if (key.compare("9") == 0)
		{
			if (value)
			{
				//Renderer->RemoveActor(inputActor);
				value = !value;
			}
			else
			{
				//Renderer->AddActor(inputActor);
				value = !value;
			}
		}

		vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		//inputMapper->SetInputData(input);
		mapper->SetInputData(updateData);
		mapper->ScalarVisibilityOff();
		//vtkSmartPointer<vtkActor> Actor =			vtkSmartPointer<vtkActor>::New();
		Actor->SetMapper(mapper);
		Actor->GetProperty()->SetColor(1, 1, 0);
		vtkSmartPointer<vtkProperty> backFaces =
			vtkSmartPointer<vtkProperty>::New();
		backFaces->SetSpecular(0.5);
		backFaces->SetDiffuse(.2);
		backFaces->SetAmbient(1.0);
		backFaces->SetAmbientColor(.4, .4, .0);
		backFaces->SetColor(.3, .3, .3);
		Actor->SetBackfaceProperty(backFaces);
		Renderer->AddActor(Actor);

		if (value)Renderer->AddActor(inputActor);
		else
			Renderer->RemoveActor(inputActor);
		Renderer->GetRenderWindow()->Render();
		Renderer->Render();
		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

	void ReallyDeletePoint(vtkSmartPointer<vtkPoints> points, vtkIdType id)
	{
		vtkSmartPointer<vtkPoints> newPoints =
			vtkSmartPointer<vtkPoints>::New();

		for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
		{
			if (i != id)
			{
				double p[3];
				points->GetPoint(i, p);
				newPoints->InsertNextPoint(p);
			}
		}

		points->ShallowCopy(newPoints);
	}
	vtkActor* Actor;
	vtkActor* inputActor;
	vtkPolyData* updateData;
	vtkPolyData* inputData;
	vtkPolyData* BooleanTargetData;
	vtkPolyData* SleeveData;
	vtkRenderer* Renderer;
};
vtkStandardNewMacro(KeyPressInteractorStyle);


int main(int, char *[])
{
	vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszPathName2 = { "D:\\target1.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();
	input = reader->GetOutput();

	vtkSmartPointer<vtkPolyData> Lineinput = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszLinePathName = { "D:\\TestPoints.vtp" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkXMLPolyDataReader> Linereader =
		vtkSmartPointer<vtkXMLPolyDataReader>::New();
	Linereader->SetFileName(lpszLinePathName);
	Linereader->Update();

	vtkSmartPointer<vtkPoints> points = Linereader->GetOutput()->GetPoints();

	vtkSmartPointer<vtkParametricSpline>mpSelectSpline = vtkSmartPointer<vtkParametricSpline>::New();
	mpSelectSpline->SetPoints(points);
	mpSelectSpline->Modified();

	vtkSmartPointer<vtkParametricFunctionSource> FunctionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
	FunctionSource->SetParametricFunction(mpSelectSpline);
	FunctionSource->Update();

	vtkSmartPointer<vtkPolyData> TmpPolydata = vtkSmartPointer<vtkPolyData>::New();

	std::cout << "Before decimation" << std::endl << "------------" << std::endl;
	std::cout << "There are " << input->GetNumberOfPoints() << " points." << std::endl;
	std::cout << "There are " << input->GetNumberOfPolys() << " polygons." << std::endl;

	vtkSmartPointer<vtkCleanPolyData> cp = vtkSmartPointer<vtkCleanPolyData>::New();
	cp->SetInputData(input);
	cp->Update();

	vtkSmartPointer<vtkPolyDataNormals> norma = vtkPolyDataNormals::New();
	norma->SetInputConnection(cp->GetOutputPort());
	norma->ComputeCellNormalsOff();
	norma->ComputePointNormalsOff();
	norma->Update();


	vtkSmartPointer<vtkSelectPolyData> mpSelectedLoopArea = vtkSmartPointer<vtkSelectPolyData>::New();
	mpSelectedLoopArea->SetInputConnection(norma->GetOutputPort());
	mpSelectedLoopArea->SetLoop(FunctionSource->GetOutput()->GetPoints()); // 루프 설정할 포인터 넘김
	mpSelectedLoopArea->GenerateSelectionScalarsOn();
	mpSelectedLoopArea->SetSelectionModeToClosestPointRegion();

	vtkSmartPointer<vtkClipPolyData> clip = vtkSmartPointer<vtkClipPolyData>::New();
	clip->SetInputConnection(mpSelectedLoopArea->GetOutputPort());
	clip->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	confilter->SetInputConnection(clip->GetOutputPort());
	confilter->SetExtractionModeToLargestRegion();
	confilter->Update();

	////// Generate normals
	vtkSmartPointer<vtkPolyDataNormals> normals =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection(confilter->GetOutputPort());
	normals->ComputeCellNormalsOn();
	normals->ComputePointNormalsOn();
	normals->Update();

	TmpPolydata->DeepCopy(normals->GetOutput());

	TmpPolydata->GetPointData()->GetNormals();

	vtkDataArray* ClipDataNormalArray = TmpPolydata->GetPointData()->GetNormals();
		
	ClipDataNormalArray = TmpPolydata->GetPointData()->GetNormals();

	std::cerr << "number of UpdataData Points :" << TmpPolydata->GetPoints()->GetNumberOfPoints() << std::endl;
	std::cerr << "number of UpdataData Normal :" << ClipDataNormalArray->GetNumberOfTuples() << std::endl;

	//for (int i = 0; i < TmpPolydata->GetPoints()->GetNumberOfPoints(); i++)
	//{
	//	//vtkP(vtkSphereSource, testSource);
	//	vtkVector3d Point(TmpPolydata->GetPoints()->GetPoint(i));
	//	vtkVector3d NormalV(ClipDataNormalArray->GetTuple(i));

	//	NormalV = NormalV * 0.4;
	//	Point = Point + NormalV;

	//	TmpPolydata->GetPoints()->SetPoint(i, Point.GetData());
	//}

	//### Create isosurface###
	vtkSmartPointer<vtkImplicitModeller> implicitModeller = vtkImplicitModeller::New();
	implicitModeller->SetSampleDimensions(50, 50, 50);
	implicitModeller->SetInputData(TmpPolydata);
	implicitModeller->AdjustBoundsOn();
	implicitModeller->SetAdjustDistance(.1); // Adjust by 10%
	implicitModeller->SetMaximumDistance(.1);

	double bounds[6];
	TmpPolydata->GetBounds(bounds);
	//implicitModeller->GetOutput()->GetBounds(bounds);
	double xrange = bounds[1] - bounds[0];

	// Create the 0 isosurface
	vtkSmartPointer<vtkContourFilter> contourFilter =
		vtkSmartPointer<vtkContourFilter>::New();
	contourFilter->SetInputConnection(implicitModeller->GetOutputPort());
	contourFilter->SetValue(0, xrange / 30.0); // 30% of xrange


	//### clip under  - 1###
	//vtkReverseSense* reverseIN = vtkReverseSense::New();;
	//reverseIN->SetInputConnection(contourFilter->GetOutputPort());
	//reverseIN->ReverseNormalsOn();
	//reverseIN->Update();
	//vtkPolyDataNormals* NormalFilter = vtkPolyDataNormals::New();
	//NormalFilter->SetInputData(reverseIN->GetOutput());
	//NormalFilter->Update();
	//
	//vtkBox* BodyClipBound = vtkBox::New();
	//double BodyBounds[6];
	//NormalFilter->GetOutput()->GetBounds(BodyBounds);
	//BodyClipBound->SetBounds(BodyBounds);

	//vtkClipPolyData* BodyBoundsClip = vtkClipPolyData::New();
	//BodyBoundsClip->SetInputData(norma->GetOutput());
	//BodyBoundsClip->SetClipFunction(BodyClipBound);
	//BodyBoundsClip->InsideOutOn();
	//BodyBoundsClip->Update();

	//vtkBooleanOperationPolyDataFilter* booleanStep1 =
	//	vtkBooleanOperationPolyDataFilter::New();
	//booleanStep1->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
	//booleanStep1->SetInputData(0, NormalFilter->GetOutput());
	//booleanStep1->SetInputData(1, BodyBoundsClip->GetOutput());
	//booleanStep1->Update();

	//### clip under - 2 ###



	//### Visualize ###
	vtkSmartPointer<vtkPolyDataMapper> lineMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	lineMapper->SetInputData(FunctionSource->GetOutput());

	vtkSmartPointer<vtkActor> lineActor =
		vtkSmartPointer<vtkActor>::New();
	lineActor->SetMapper(lineMapper);
	lineActor->GetProperty()->SetColor(1., 0., 0.);

	vtkSmartPointer<vtkPolyDataMapper> inputMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	//inputMapper->SetInputData(TmpPolydata);
	inputMapper->SetInputConnection(contourFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);

	vtkSmartPointer<vtkPolyDataMapper> decimatedMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	decimatedMapper->SetInputConnection(reader->GetOutputPort());
	vtkSmartPointer<vtkActor> decimatedActor =
		vtkSmartPointer<vtkActor>::New();
	decimatedActor->SetMapper(decimatedMapper);
	decimatedActor->GetProperty()->SetColor(1, 1, .5);

	// There will be one render window
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(1800, 900);

	// Create a camera for all renderers
	vtkSmartPointer<vtkCamera> camera =
		vtkSmartPointer<vtkCamera>::New();

	// Define viewport ranges
	// (xmin, ymin, xmax, ymax)
	double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
	double rightViewport[4] = { 0.5, 0.0, 1.0, 1.0 };

	// Setup both renderers
	vtkSmartPointer<vtkRenderer> leftRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(leftRenderer);
	leftRenderer->SetViewport(leftViewport);
	leftRenderer->SetBackground(.6, .5, .4);
	leftRenderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderer> rightRenderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(rightRenderer);
	rightRenderer->SetViewport(rightViewport);
	rightRenderer->SetBackground(.4, .5, .6);
	rightRenderer->SetActiveCamera(camera);

	// Add the sphere to the left and the cube to the right
	leftRenderer->AddActor(lineActor);
	leftRenderer->AddActor(inputActor);
	
	rightRenderer->AddActor(lineActor);
	rightRenderer->AddActor(decimatedActor);
	rightRenderer->AddActor(inputActor);
	//rightRenderer->AddActor(warpactor);

	leftRenderer->ResetCamera();

	rightRenderer->ResetCamera();
	// And one interactor
	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	vtkSmartPointer<KeyPressInteractorStyle> style =
		vtkSmartPointer<KeyPressInteractorStyle>::New();
	style->Actor = decimatedActor;
	style->inputActor = inputActor;
	style->inputData = reader->GetOutput();
	vtkSmartPointer<vtkPolyData> update = vtkSmartPointer<vtkPolyData>::New();
	update->DeepCopy(reader->GetOutput());
	style->updateData = update;
	style->Renderer = rightRenderer;
	interactor->SetInteractorStyle(style);
	style->SetCurrentRenderer(rightRenderer);

	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}