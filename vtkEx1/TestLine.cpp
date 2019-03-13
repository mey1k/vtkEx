#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFillHolesFilter.h>
#include <vtkTransform.h>
#include <vtkCutter.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPlane.h>
#include <thread>
#include <chrono>
#include <vtkTransformPolyDataFilter.h>
#include <vtkBox.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPointLocator.h>
#include <vtkCellPicker.h>
#include <vtkKdTree.h>
#include <vtkKdTreePointLocator.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkQuadricClustering.h>
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
#include <vtkRenderWindow.h>
#include <vtkVector.h>
#include <vtkVectorOperators.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSortDataArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkWarpVector.h>
#include <vtkDataSetAttributes.h>
#include <vtkCleanPolyData.h>
#include <vtkStripper.h>
#include <vtkActor.h>
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
#include <vtkCellLocator.h>
#include <vtkLineSource.h>
#include <vtkPolyLine.h>
#include <vtkOBBTree.h>
#include <vtkMath.h>
#include <vtkModifiedBSPTree.h>
#include <vtkNew.h>
//#include <vtkCollisionDetectionFilter>

#include <vtkVoxelContoursToSurfaceFilter.h>
#include <vtkStripper.h>
#include <vtkTransformFilter.h>
#include <vtkFloatArray.h>
#include <vtkWeightedTransformFilter.h>
#include <vtkArrowSource.h>
#include <vtkOutlineFilter.h>
#include <vtkWarpTo.h>
#include <vtkDelaunay2D.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetMapper.h>
#include <vtkSelectPolyData.h>
#include <vtkLODActor.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkExtractSelection.h>
#include <vtkInformation.h> // property->set
#include <vtkExtractEdges.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetSurfaceFilter.h>

bool GetPointNormals(vtkPolyData* polydata);

vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> input2 = vtkSmartPointer<vtkPolyData>::New();

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);

//vtkSmartPointer<vtkActor> contoursActor = vtkSmartPointer<vtkActor>::New();
//vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	bool value = 1;
	static KeyPressInteractorStyle* New();

	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

	void AddRenderer(vtkPolyData* polydata)
	{
		vtkSmartPointer<vtkPolyDataMapper> Mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
		Mapper->SetInputData(polydata);
		//Mapper->SetResolveCoincidentTopologyToPolygonOffset();
		Mapper->SetResolveCoincidentTopologyToShiftZBuffer();
		//Mapper->ScalarVisibilityOff();
		vtkSmartPointer<vtkActor> Actor =
			vtkSmartPointer<vtkActor>::New();
		Actor->SetMapper(Mapper);
		Actor->GetProperty()->SetColor(1, 0, 0);
		Actor->GetProperty()->SetLineWidth(2.0);

		Renderer->AddActor(outActor);
	}
	virtual void OnKeyPress()
	{
		// Get the keypress
		std::string key = this->Interactor->GetKeySym();
		Renderer->RemoveActor(Actor);
		std::cerr << " - ";

		if (key.compare("0") == 0)
		{
			updateData->DeepCopy(inputData);
			std::cerr << "Do Initialize" << std::endl;
		}// "s" for "s"elect
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
		else if (key.compare("2") == 0)
		{
			///알고리즘 전제조건 : 단일 라인으로 이루어진 Contour가 검출되어야함
			vtkFeatureEdges* featureEdges = vtkFeatureEdges::New();
			//featureEdges->SetInputConnection(updateData->GetOutputPort());
			featureEdges->SetInputData(updateData);
			featureEdges->BoundaryEdgesOff();
			featureEdges->FeatureEdgesOn();
			featureEdges->ManifoldEdgesOff();
			featureEdges->NonManifoldEdgesOff();
			featureEdges->Update();
			AddRenderer(featureEdges->GetOutput());
		}
		else if (key.compare("4") == 0)
		{
			vtkSmartPointer<vtkPlane> plane =
				vtkSmartPointer<vtkPlane>::New();
			plane->SetOrigin(updateData->GetCenter());
			plane->SetNormal(-0.5, 1, 0);

			std::cout << "updateDataPoints" << updateData->GetNumberOfPoints() << std::endl;

			double minBound[3];
			minBound[0] = updateData->GetBounds()[0];
			minBound[1] = updateData->GetBounds()[2];
			minBound[2] = updateData->GetBounds()[4];

			double maxBound[3];
			maxBound[0] = updateData->GetBounds()[1];
			maxBound[1] = updateData->GetBounds()[3];
			maxBound[2] = updateData->GetBounds()[5];

			double center[3];
			center[0] = updateData->GetCenter()[0];
			center[1] = updateData->GetCenter()[1];
			center[2] = updateData->GetCenter()[2];

			double distanceMin = sqrt(vtkMath::Distance2BetweenPoints(minBound, center));
			double distanceMax = sqrt(vtkMath::Distance2BetweenPoints(maxBound, center));

			cout << minBound[0] << " " << maxBound[0];
			cout << minBound[1] << " " << maxBound[1];
			cout << minBound[2] << " " << maxBound[2];


			vtkSmartPointer<vtkPoints> undercutPoint =
				vtkSmartPointer<vtkPoints>::New();

			//// Create cutter
			//vtkSmartPointer<vtkCutter> cutter =
			//	vtkSmartPointer<vtkCutter>::New();
			//cutter->SetCutFunction(plane);
			//cutter->SetInputData(updateData);
			//cutter->GenerateValues(200, -distanceMin, distanceMax);
			//cutter->Update();

			//vtkSmartPointer<vtkCutter> copyCutter = vtkSmartPointer<vtkCutter>::New();
			//copyCutter->SetCutFunction(plane);
			//copyCutter->SetInputData(updateData);
			//cutter->GenerateValues(200, -distanceMin, distanceMax);
			//copyCutter->Update();

			//vtkSmartPointer<vtkPolyDataNormals> normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
			//normalFilter->SetInputData(cutter->GetOutput());
			//normalFilter->ComputePointNormalsOn();
			//normalFilter->Update();
			//normalFilter->GetOutput()->GetPoint((vtkIdType)1);

			vtkSmartPointer<vtkPoints> meshPoints = vtkSmartPointer<vtkPoints>::New();
			meshPoints = updateData->GetPoints();
			
			cout << "meshPointsCount" << meshPoints->GetNumberOfPoints() << std::endl;

			//vtkIdType NowPointID, PrePointID;
			//vtkSmartPointer<vtkIdList> cellIdList;
			//vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

			//for (NowPointID = 0; NowPointID < cutter->GetOutput()->GetPoints()->GetNumberOfPoints(); NowPointID++)
			//{
			//	cellIdList =	vtkSmartPointer<vtkIdList>::New();
			//	cutter->GetOutput()->GetPointCells(NowPointID, cellIdList);

			//	if (cellIdList->GetNumberOfIds() == 1)
			//	{
			//		//std::cerr << "@NowPointID :" << NowPointID <<" | Start or end"<< std::endl;
			//		std::cerr << "@cellID:" << cellIdList ->GetId(0)<<" | Start or end CELL"<< std::endl;
			//		connectedVertices->InsertNextId(NowPointID);
			//		PrePointID = NowPointID;
			//		break;
			//	}
			//}

			/*vtkSmartPointer<vtkPoints> ReMakePoints = vtkSmartPointer<vtkPoints>::New();
			vtkSmartPointer<vtkCellArray> ReMakeCells = vtkSmartPointer<vtkCellArray>::New();
			int numCell = 0;
			ReMakePoints->InsertNextPoint(cutter->GetOutput()->GetPoints()->GetPoint(NowPointID));*/

			// Create the locator
			/*vtkSmartPointer<vtkCellLocator> mpPointLocator =
				vtkSmartPointer<vtkCellLocator>::New();
			mpPointLocator->SetDataSet(cutter->GetOutput());
			mpPointLocator->BuildLocator();*/

			/*int cutterNumberOfPoints = cutter->GetOutput()->GetPoints()->GetNumberOfPoints();*/

			double tolerance = 0.001;
			double pos1[] = { 0.0, 0.0, 0.0 };
			double pos2[] = { 0.0, 0.0, 0.0 };

			double prepos1[] = { 0.0, 0.0, 0.0 };
			double prepos2[] = { 0.0, 0.0, 0.0 };

			double pcoords[] = { 0.0, 0.0, 0.0 };
			int subId = 0;
			double t = 0;
			int CurrentPointCnt = 1;
			
			vtkVector3d cutterDirection, startPoint, endPoint;
			//start
			//cutter->GetOutput()->GetPoints()->GetPoint(0, startPoint.GetData());
			////end
			//cutter->GetOutput()->GetPoints()->GetPoint(cutterNumberOfPoints-1, endPoint.GetData());
			//direction
			vtkMath::Subtract(endPoint.GetData(), startPoint.GetData(), cutterDirection.GetData());
			vtkMath::Normalize(cutterDirection.GetData());

			vtkVector3d insersionVec;
			insersionVec[0] = 0;
			insersionVec[1] = 0;
			insersionVec[2] = maxBound[2];

			vtkVector3d cutterPoint;
			vtkVector3d topPoint, bottomPoint,underPoint;

			bool isFalt = false;

			double dValue = 1.05;
			double value = floor(100.*(dValue + 0.005)) / 100.;
			int test = 0;

			int insertionCount = 0;

			vtkSmartPointer<vtkIdTypeArray> ids =
				vtkSmartPointer<vtkIdTypeArray>::New();
			ids->SetNumberOfComponents(1);

			vtkSmartPointer<vtkIdTypeArray> undercutids =
				vtkSmartPointer<vtkIdTypeArray>::New();
			undercutids->SetNumberOfComponents(1);

			/*bool hasPointNormals = GetPointNormals(updateData);

			std::cout << "haspointnormals?..." << hasPointNormals << std::endl;*/
			cout << "updateData1" << updateData->GetNumberOfPoints() << std::endl;

			vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION <= 5
			normalGenerator->SetInput(polydata);
#else
			normalGenerator->SetInputData(updateData);
#endif
			normalGenerator->ComputePointNormalsOn();
			normalGenerator->ComputeCellNormalsOff();
			normalGenerator->Update();

			updateData = normalGenerator->GetOutput();

			vtkSmartPointer<vtkPolyData> originalData = vtkSmartPointer<vtkPolyData>::New();

			originalData->DeepCopy(updateData);

			vtkSmartPointer<vtkCellLocator> mpPointLocator =
				vtkSmartPointer<vtkCellLocator>::New();
			mpPointLocator->SetDataSet(updateData);
			mpPointLocator->BuildLocator();

			bool isok = false;

			int wrong = 0;
			int wrongtest = 0;

			double normalizedZ[3];
			normalizedZ[0] = insersionVec[0];
			normalizedZ[1] = insersionVec[1];
			normalizedZ[2] = insersionVec[2];

			vtkMath::Normalize(normalizedZ);

			for (vtkIdType i = 0; i < updateData->GetNumberOfPoints(); i++)
			{
				double p[3];

				//while (!isok)
				//{
					updateData->GetPoint(i, p);

					topPoint = (vtkVector3d)p + (vtkVector3d)normalizedZ*10;
					//bottomPoint = (vtkVector3d)p - insersionVec;
					bottomPoint = (vtkVector3d)p;

					double point0[3];
					double point1[3];

					point0[0] = topPoint[0];
					point0[1] = topPoint[1];
					point0[2] = topPoint[2]; 

					point1[0] = bottomPoint[0];
					point1[1] = bottomPoint[1];
					point1[2] = bottomPoint[2];

					mpPointLocator->IntersectWithLine(point0, point1, tolerance, t, pos1, pcoords, subId);
					mpPointLocator->IntersectWithLine(point1, point0, tolerance, t, pos2, pcoords, subId);

					//if (CompareDoubleAbsoulte(pos1, prepos1) && CompareDoubleAbsoulte(pos2, prepos2)) break;

					if (!CompareDoubleAbsoulte(pos1, pos2))
					{
						undercutids->InsertNextValue(i);
					
						undercutPoint->InsertPoint(i, p);

						/*double normalizedX[3];
						double normalizedY[3];
						double testPoint[3] = { p[0],p[1],p[2] };

						vtkDataArray* normalsGeneric = updateData->GetPointData()->GetNormals();
						normalsGeneric->GetTuple(i, testPoint);

						double arbitrary[3];
						arbitrary[0] = testPoint[0];
						arbitrary[1] = testPoint[1];
						arbitrary[2] = testPoint[2];
						vtkMath::Normalize(arbitrary);

						testPoint[0] += arbitrary[0];
						testPoint[1] += arbitrary[1];
						testPoint[2] += arbitrary[2];

						vtkMath::Cross(testPoint, normalizedZ, normalizedX);
						vtkMath::Normalize(normalizedX);

						vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
						vtkMath::Normalize(normalizedY);

						vtkVector3d test;

						test = (vtkVector3d)p + (vtkVector3d)normalizedY*0.5;

						updateData->GetPoints()->SetPoint(i, test.GetData());*/

						/*double testPoint[3] = { p[0],p[1],p[2] };
						vtkVector3d test;

						vtkDataArray* normalsGeneric = updateData->GetPointData()->GetNormals();
						normalsGeneric->GetTuple(i, testPoint);

						vtkMath::Normalize(testPoint);

						test = (vtkVector3d)p + (vtkVector3d)testPoint*0.5;

						updateData->GetPoints()->SetPoint(i, test.GetData());*/

						//wrong++;

						//if (wrong == 1000) break;

						//ReallyDeletePoint(updateData->GetPoints(), i);

						/*updateData->RemoveDeletedCells();
						updateData->DeleteLinks();
						updateData->Modified();*/

						/*if (wrongtest % 100 == 0)
						{

							vtkSmartPointer<vtkLineSource> lineSource =
								vtkSmartPointer<vtkLineSource>::New();
							lineSource->SetPoint1(point0);
							lineSource->SetPoint2(point1);
							lineSource->Update();

							vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
							mapper->SetInputConnection(lineSource->GetOutputPort());

							vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
							actor->SetMapper(mapper);
							actor->GetProperty()->SetColor(colorRGB(0, 0, 255));
							actor->GetProperty()->SetOpacity(0.7);

							Renderer->AddActor(actor);
						}
						wrongtest++;
						wrong++;

						if (wrong == 50) break;*/

						//vtkSTLWriter* stlWriter = vtkSTLWriter::New();
						//const char* lpszPathName = { "D:\\falut.stl" };
						//stlWriter->SetFileName(lpszPathName);
						//stlWriter->SetInputData(updateData);
						////stlWriter->SetInputConnection(updateData->GetOutputPort());
						//stlWriter->Write();
						////wrong++;

						vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(updateData, i);

						for (int j = 0; j < connectedVertices->GetNumberOfIds(); j++)
						{
							/*updateData->GetPoint(connectedVertices->GetId(j), p);

							double normalizedX[3];
							double normalizedY[3];
							double testPoint[3] = { p[0],p[1],p[2] };

							vtkDataArray* normalsGeneric = updateData->GetPointData()->GetNormals();
							normalsGeneric->GetTuple(connectedVertices->GetId(j), testPoint);

							double arbitrary[3];
							arbitrary[0] = testPoint[0];
							arbitrary[1] = testPoint[1];
							arbitrary[2] = testPoint[2];
							vtkMath::Normalize(arbitrary);

							testPoint[0] += arbitrary[0];
							testPoint[1] += arbitrary[1];
							testPoint[2] += arbitrary[2];

							vtkMath::Cross(testPoint, normalizedZ, normalizedX);
							vtkMath::Normalize(normalizedX);

							vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
							vtkMath::Normalize(normalizedY);

							vtkVector3d test;

							test = (vtkVector3d)p + (vtkVector3d)normalizedY*0.5;

							updateData->GetPoints()->SetPoint(connectedVertices->GetId(j), test.GetData());*/

							undercutids->InsertNextValue(connectedVertices->GetId(j));
						}
					}
					else
					{
						wrong = 0;
						isok = true;
						ids->InsertNextValue(i);
					}
				//}

				wrong = 0;
				isok = false;
			}

			/*vtkSmartPointer<vtkSelectPolyData> loop =
				vtkSmartPointer<vtkSelectPolyData>::New();
			loop->SetInputData(updateData);
			loop->SetLoop(undercutPoint);
			loop->GenerateSelectionScalarsOn();
			loop->SetSelectionModeToSmallestRegion();

			vtkSmartPointer<vtkClipPolyData> clip =
				vtkSmartPointer<vtkClipPolyData>::New();
			clip->SetInputConnection(loop->GetOutputPort());*/		

			//vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			//const char* lpszPathName = { "D:\\Result1.stl" };
			//stlWriter->SetFileName(lpszPathName);
			////stlWriter->SetInputData(updateData);
			//stlWriter->SetInputConnection(cleaner->GetOutputPort());
			//stlWriter->Write();
			//vtkSmartPointer<vtkPolyDataMapper> updateMapper =
			//	vtkSmartPointer<vtkPolyDataMapper>::New();
			////updateMapper->SetInputData(updateData);
			//updateMapper->SetInputConnection(clip->GetOutputPort());

			//vtkSmartPointer<vtkActor> updateActor =
			//	vtkSmartPointer<vtkActor>::New();
			//updateActor -> SetMapper(updateMapper);
			//updateActor->GetProperty()->SetColor(colorRGB(255, 0, 255));

			//Renderer->AddActor(updateActor);

			vtkSmartPointer<vtkSelectionNode> selectionNode =
				vtkSmartPointer<vtkSelectionNode>::New();
			selectionNode->SetFieldType(vtkSelectionNode::POINT);
			selectionNode->SetContentType(vtkSelectionNode::INDICES);
			selectionNode->SetSelectionList(undercutids);
			selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

			vtkSmartPointer<vtkSelection> selection =
				vtkSmartPointer<vtkSelection>::New();
			selection->AddNode(selectionNode);

			vtkSmartPointer<vtkExtractSelection> extractSelection =
				vtkSmartPointer<vtkExtractSelection>::New();

			extractSelection->SetInputData(0, updateData);
#if VTK_MAJOR_VERSION <= 5
			extractSelection->SetInput(1, selection);
#else
			extractSelection->SetInputData(1, selection);
#endif
			extractSelection->Update();



			vtkSmartPointer<vtkUnstructuredGrid> notSelected =
				vtkSmartPointer<vtkUnstructuredGrid>::New();
			notSelected->ShallowCopy(extractSelection->GetOutput());

			vtkSmartPointer<vtkExtractEdges> extractEdges =
				vtkSmartPointer<vtkExtractEdges>::New();
			//extractEdges->SetInputConnection(notSelected->GetOutputPort());
			extractEdges->SetInputData(notSelected);
			extractEdges->Update();

			vtkSmartPointer<vtkPolyDataMapper> ex =
				vtkSmartPointer<vtkPolyDataMapper> ::New();
			ex->SetInputConnection(extractEdges->GetOutputPort());

			vtkSmartPointer<vtkActor> exac =
				vtkSmartPointer<vtkActor>::New();
			exac->SetMapper(ex);
			exac->GetProperty()->SetColor(colorRGB(255, 0, 255));

			vtkSmartPointer<vtkDataSetMapper> notSelectedMapper =
				vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
			notSelectedMapper->SetInputConnection(notSelected->GetProducerPort());
#else
			notSelectedMapper->SetInputData(notSelected);
#endif
			vtkSmartPointer<vtkActor> notSelectedActor =
				vtkSmartPointer<vtkActor>::New();
			notSelectedActor->SetMapper(notSelectedMapper);
			notSelectedActor->GetProperty()->SetColor(colorRGB(0, 255, 255));

			vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(extractEdges->GetOutputPort());
			vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(colorRGB(0, 255, 255));

			Renderer->AddActor(notSelectedActor);


			vtkSmartPointer<vtkDataSetSurfaceFilter> dsf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
			dsf->SetInputData(notSelected);
			dsf->Update();

			/*vtkSmartPointer<vtkTriangleFilter> triangleFilterx = vtkSmartPointer<vtkTriangleFilter>::New();
			triangleFilterx->
			triangleFilterx->SetInputConnection(dsf->GetOutputPort());
			triangleFilterx->Update();*/

//			vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
//				vtkSmartPointer<vtkFillHolesFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//			fillHolesFilter->SetInputConnection(input->GetProducerPort());
//#else
//			fillHolesFilter->SetInputConnection(dsf->GetOutputPort());
//#endif
//			fillHolesFilter->SetHoleSize(100.0);
//
//			// Make the triangle windong order consistent
//			vtkSmartPointer<vtkPolyDataNormals> normals =
//				vtkSmartPointer<vtkPolyDataNormals>::New();
//			normals->SetInputConnection(fillHolesFilter->GetOutputPort());
//			normals->ConsistencyOn();
//			normals->SplittingOff();
//			normals->Update();
//
//			// Restore the original normals
//			normals->GetOutput()->GetPointData()->
//				SetNormals(dsf->GetOutput()->GetPointData()->GetNormals());
//
// Create the outline
			vtkSmartPointer<vtkOutlineFilter> outline =
				vtkSmartPointer<vtkOutlineFilter>::New();
#if VTK_MAJOR_VERSION <= 5
			outline->SetInput(sphere);
#else
			outline->SetInputData(dsf->GetOutput());
			//outline->SetInputConnection(dsf->GetOutput());
#endif
			vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			outlineMapper->SetInputConnection(outline->GetOutputPort());
			vtkSmartPointer<vtkActor> outlineActor =
				vtkSmartPointer<vtkActor>::New();
			outlineActor->SetMapper(outlineMapper);
			outlineActor->GetProperty()->SetColor(0, 0, 0);


			vtkSTLWriter* stlWriter = vtkSTLWriter::New();
			const char* lpszPathName = { "D:\\Result1.stl" };
			stlWriter->SetFileName(lpszPathName);
			//stlWriter->SetInputData(notSelected);
			stlWriter->SetInputConnection(dsf->GetOutputPort());
			stlWriter->Write();

			vtkSTLWriter* stlWriter2 = vtkSTLWriter::New();
			const char* lpszPathName2 = { "D:\\Result2.stl" };
			stlWriter2->SetFileName(lpszPathName2);
			stlWriter2->SetInputData(updateData);
			stlWriter2->Write();


			//Renderer->AddActor(exac);
			
			//Renderer->AddActor(outlineActor);

			//vtkSmartPointer<vtkSelectPolyData> loop =
			//	vtkSmartPointer<vtkSelectPolyData>::New();
			//loop->SetInputData(updateData);
			////loop->SetInputConnection(updateData->GetOutputPort());
			//loop->SetLoop(selectionPoints);
			//loop->GenerateSelectionScalarsOn();
			//loop->SetSelectionModeToSmallestRegion(); //negative scalars inside

			//vtkSmartPointer<vtkClipPolyData> clip = //clips out positive region
			//	vtkSmartPointer<vtkClipPolyData>::New();
			//clip->SetInputConnection(loop->GetOutputPort());

			//vtkSmartPointer<vtkPolyDataMapper> clipMapper =
			//	vtkSmartPointer<vtkPolyDataMapper>::New();
			//clipMapper->SetInputConnection(clip->GetOutputPort());

			//vtkSmartPointer<vtkLODActor> clipActor =
			//	vtkSmartPointer<vtkLODActor>::New();
			//clipActor->SetMapper(clipMapper);

			//Renderer->AddActor(clipActor);

			// Add the grid points to a polydata object
			/*vtkSmartPointer<vtkPolyData> polydata =
				vtkSmartPointer<vtkPolyData>::New();
			polydata->SetPoints(meshPoints);

			vtkSmartPointer<vtkCleanPolyData> cleaner =
				vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->SetInputData(polydata);

			vtkSmartPointer<vtkPolyDataMapper> delaunayMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			delaunayMapper->SetInputConnection(cleaner->GetOutputPort());

			vtkSmartPointer<vtkActor> delaunayActor =
				vtkSmartPointer<vtkActor>::New();
			delaunayActor->SetMapper(delaunayMapper);
			delaunayActor->GetProperty()->SetColor(1, 1, 1);*/
			//cleaner->SetInputConnection(polydata->GetOutputPort());

			// Generate a tetrahedral mesh from the input points. By
			// default, the generated volume is the convex hull of the points.
			/*vtkSmartPointer<vtkDelaunay3D> delaunay3D =
				vtkSmartPointer<vtkDelaunay3D>::New();
			delaunay3D->SetInputConnection(cleaner->GetOutputPort());

			vtkSmartPointer<vtkDataSetMapper> delaunayMapper =
				vtkSmartPointer<vtkDataSetMapper>::New();
			delaunayMapper->SetInputConnection(delaunay3D->GetOutputPort());

			vtkSmartPointer<vtkActor> delaunayActor =
				vtkSmartPointer<vtkActor>::New();
			delaunayActor->SetMapper(delaunayMapper);
			delaunayActor->GetProperty()->SetColor(1, 0, 0);*/

			//Renderer->AddActor(delaunayActor);

			/*vtkSmartPointer<vtkPolyData> polydata =
				vtkSmartPointer<vtkPolyData>::New();
			polydata->SetPoints(meshPoints);

			vtkSmartPointer<vtkDelaunay3D> delaunay3D =
				vtkSmartPointer<vtkDelaunay3D>::New();
			delaunay3D->SetInputData(polydata);

			vtkSmartPointer<vtkSurfaceReconstructionFilter> surf =
				vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();*/

			//surf->SetInput(polydata);

			//surf->SetInputData(polydata);

			/*vtkSmartPointer<vtkCleanPolyData> cleaner =
				vtkSmartPointer<vtkCleanPolyData>::New();
			cleaner->SetInputData(polydata);*/

//			// Triangulate the grid points
//			vtkSmartPointer<vtkDelaunay2D> delaunay =
//				vtkSmartPointer<vtkDelaunay2D>::New();
//#if VTK_MAJOR_VERSION <= 5
//			delaunay->SetInput(polydata);
//#else
//			delaunay->SetInputData(polydata);
//#endif
//			delaunay->Update();
//
//			vtkSmartPointer<vtkWarpTo> warpTo =
//				vtkSmartPointer<vtkWarpTo>::New();
//			warpTo->SetInputConnection(delaunay->GetOutputPort());
//			warpTo->SetPosition(0, 0, 0);
//			warpTo->SetScaleFactor(0);
//			warpTo->AbsoluteOn();

			//vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
			//	vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
			//smoothFilter->SetInputConnection(delaunay->GetOutputPort());
			//smoothFilter->SetNumberOfIterations(15);
			//smoothFilter->SetRelaxationFactor(0.1);
			//smoothFilter->FeatureEdgeSmoothingOff();
			//smoothFilter->BoundarySmoothingOn();
			//smoothFilter->Update();

			//// Update normals on newly smoothed polydata
			//vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
			//normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
			//normalGenerator->ComputePointNormalsOn();
			//normalGenerator->ComputeCellNormalsOn();
			//normalGenerator->Update();

			/*vtkSmartPointer<vtkContourFilter> cf =
				vtkSmartPointer<vtkContourFilter>::New();
			cf->SetInputConnection(surf->GetOutputPort());
			cf->SetValue(0, 0.0);

			vtkSmartPointer<vtkPolyDataMapper> inputMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			inputMapper->SetInputConnection(delaunay3D->GetOutputPort());
			vtkSmartPointer<vtkActor> inputActor =
				vtkSmartPointer<vtkActor>::New();
			inputActor->SetMapper(inputMapper);
			inputActor->GetProperty()->SetColor(colorRGB(0, 255, 255));*/


			/*vtkSmartPointer<vtkPolyDataMapper> smoothedMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			smoothedMapper->SetInputConnection(normalGenerator->GetOutputPort());
			vtkSmartPointer<vtkActor> smoothedActor =
				vtkSmartPointer<vtkActor>::New();
			smoothedActor->SetMapper(smoothedMapper);
			smoothedActor->GetProperty()->SetColor(colorRGB(0, 255, 255));*/

			//Renderer->AddActor(inputActor);
//
//			// There will be one render window
//			vtkSmartPointer<vtkRenderWindow> renderWindow =
//				vtkSmartPointer<vtkRenderWindow>::New();
//			renderWindow->SetSize(600, 300);
//
//			// And one interactor
//			vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//				vtkSmartPointer<vtkRenderWindowInteractor>::New();
//			interactor->SetRenderWindow(renderWindow);
//
//			// Define viewport ranges
//			// (xmin, ymin, xmax, ymax)
//			double leftViewport[4] = { 0.0, 0.0, 0.5, 1.0 };
//			double rightViewport[4] = { 0.5, 0.0, 1.0, 1.0 };
//
//			// Setup both renderers
//			vtkSmartPointer<vtkRenderer> leftRenderer =
//				vtkSmartPointer<vtkRenderer>::New();
//			renderWindow->AddRenderer(leftRenderer);
//			leftRenderer->SetViewport(leftViewport);
//			leftRenderer->SetBackground(.6, .5, .4);
//
//			// Add the input parabola to the left and the smoothed parabola to the right
//			leftRenderer->AddActor(delaunay);
//
//			renderWindow->Render();
//

			cout << "insertioncount : " << insertionCount;

//			//contourToSurface Test
//
//			vtkSmartPointer<vtkStripper> cutStrips = vtkSmartPointer<vtkStripper>::New();
//			cutStrips->SetInputConnection(cutter->GetOutputPort());
//			cutStrips->Update();
//
//			vtkSmartPointer<vtkPolyData> cutPoly = vtkSmartPointer<vtkPolyData>::New();
//			cutPoly->SetPoints((cutStrips->GetOutput())->GetPoints());
//			cutPoly->SetPolys((cutStrips->GetOutput())->GetLines());
//
//			vtkPolyData* contours = cutPoly;
//			contours->GetBounds(bounds);
//			double origin[3] = { bounds[0], bounds[2], bounds[4] };
//			double spacing[3] = { (bounds[1] - bounds[0]) /200 ,
//				(bounds[3] - bounds[2]) /200,
//				(bounds[5] - bounds[4]) };
//
//			vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
//			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//
//			vtkPoints* contourPoints = contours->GetPoints();
//			int numPoints = contourPoints->GetNumberOfPoints();
//			points->SetNumberOfPoints(numPoints);
//			for (int i = 0; i < numPoints; ++i)
//			{
//				double pt[3];
//				contourPoints->GetPoint(i, pt);
//				pt[0] = static_cast<int>((pt[0] - origin[0]) / spacing[0] + 0.5);
//				pt[1] = static_cast<int>((pt[1] - origin[1]) / spacing[1] + 0.5);
//				pt[2] = static_cast<int>((pt[2] - origin[2]) / spacing[2] + 0.5);
//				points->SetPoint(i, pt);
//			}
//			poly->SetPolys(contours->GetPolys());
//			poly->SetPoints(points);
//
//			
//
//			// Create the contour to surface filter
//			//
//			vtkSmartPointer<vtkVoxelContoursToSurfaceFilter> contoursToSurface = vtkSmartPointer<vtkVoxelContoursToSurfaceFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//			contoursToSurface->SetInput(poly);
//#else
//			contoursToSurface->SetInputData(poly);
//#endif
//			contoursToSurface->SetSpacing(spacing[0], spacing[1], spacing[2]);
//			contoursToSurface->Update();
//
//			// Rescale the output back into world coordinates and center it
//			//
//			double scaleCenter[3];
//			contoursToSurface->GetOutput()->GetCenter(scaleCenter);
//			double scaleBounds[6];
//			contoursToSurface->GetOutput()->GetBounds(scaleBounds);
//			//double center[3];
//			contours->GetCenter(center);
//
//			vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//			transformFilter->SetInputConnection(contoursToSurface->GetOutputPort());
//			
//			
//			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
//			transformFilter->SetTransform(transform);
//			transform->Translate(-scaleCenter[0], -scaleCenter[1], -scaleCenter[2]);
//			transform->Scale(
//				(bounds[1] - bounds[0]) / (scaleBounds[1] - scaleBounds[0]),
//				(bounds[3] - bounds[2]) / (scaleBounds[3] - scaleBounds[2]),
//				(bounds[5] - bounds[4]) / (scaleBounds[5] - scaleBounds[4]));
//			transform->Translate(center[0], center[1], center[2]);
//
//			// Visualize the contours
//			//
//			vtkSmartPointer<vtkPolyDataMapper> contoursMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//#if VTK_MAJOR_VERSION <= 5
//			contoursMapper->SetInput(contours);
//#else
//			contoursMapper->SetInputData(contours);
//#endif
//			contoursMapper->ScalarVisibilityOff();
//
//			vtkSmartPointer<vtkActor> contoursActor = vtkSmartPointer<vtkActor>::New();
//			contoursActor->SetMapper(contoursMapper);
//			contoursActor->GetProperty()->SetRepresentationToWireframe();
//			contoursActor->GetProperty()->ShadingOff();
//
//			// Visualize the surface
//			//
//			vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//			surfaceMapper->SetInputConnection(transformFilter->GetOutputPort());
//			surfaceMapper->ScalarVisibilityOff();
//			surfaceMapper->ImmediateModeRenderingOn();
//
//			vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
//			surfaceActor->SetMapper(surfaceMapper);
//			surfaceActor->GetProperty()->SetRepresentationToWireframe();
//			surfaceActor->GetProperty()->ShadingOff();
//
//			// Create two renderers side by side to show the contours and the surface separately
//			//
//			// Press 't' for trackball interaction
//			// Press 'r' to reset the camera
//			// Press 'w' for wireframe representation
//			// Press 's' for surface representation
//			//
//			vtkSmartPointer<vtkRenderer> renderer1 = vtkSmartPointer<vtkRenderer>::New();
//			renderer1->SetViewport(0., 0., 0.5, 1.);
//			renderer1->SetBackground(0.2, 0.2, 0.8);
//
//			vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();
//			renderer2->SetViewport(0.5, 0., 1., 1.);
//			renderer2->SetBackground(0.8, 0.2, 0.2);
//
//			vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//			renderWindow->SetSize(800, 400);
//
//			renderWindow->AddRenderer(renderer1);
//			renderWindow->AddRenderer(renderer2);
//
//			vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//			interactor->SetRenderWindow(renderWindow);
//
//			renderer1->AddViewProp(contoursActor);
//			renderer2->AddViewProp(surfaceActor);
//			renderWindow->Render();
//
//			interactor->Start();

			
			//while (1)
			//{

	

			//	copyCutter->GetOutput()->GetPoints()->GetPoint(CurrentPointCnt - 1, cutterPoint.GetData());

			//	topPoint = cutterPoint + insersionVec;
			//	bottomPoint = cutterPoint - insersionVec;

			//	double point0[3];
			//	double point1[3];

			//	point0[0] = topPoint[0];
			//	point0[1] = topPoint[1];
			//	point0[2] = topPoint[2];

			//	point1[0] = bottomPoint[0];
			//	point1[1] = bottomPoint[1];
			//	point1[2] = bottomPoint[2];

			//	mpPointLocator->IntersectWithLine(point0, point1, tolerance, t, pos1, pcoords, subId);
			//	mpPointLocator->IntersectWithLine(point1, point0, tolerance, t, pos2, pcoords, subId);

			//	//new Line
			//	if (!CompareDoubleAbsoulte(pos1, pos2))
			//	{
			//		isFalt = true;

			//		int wrongPointNum = cutter->GetOutput()->FindPoint(pos2[0], pos2[1], pos2[2]);
			//		
			//		
			//		vtkVector3d wrongPointVec = (vtkVector3d)pos2;
			//		vtkVector3d modifiedVec;

			//		if (wrongPointNum < cutterNumberOfPoints / 2)
			//		{
			//			modifiedVec = (wrongPointVec + cutterDirection) - cutterDirection * value;
			//		}
			//		else
			//		{
			//			modifiedVec = (wrongPointVec - cutterDirection) + cutterDirection * value;
			//		}

			//		pos2[0] = modifiedVec[0];
			//		pos2[1] = modifiedVec[1];
			//		pos2[2] = modifiedVec[2];

			//		cutter->GetOutput()->GetPoints()->SetPoint(wrongPointNum, pos2);
			//		

			//		/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			//		sphereSource->SetCenter(pos2);
			//		sphereSource->SetRadius(0.1);

			//		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			//		mapper->SetInputConnection(sphereSource->GetOutputPort());

			//		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			//		actor->SetMapper(mapper);
			//		actor->GetProperty()->SetColor(colorRGB(0, 0, 255));
			//		actor->GetProperty()->SetOpacity(0.7);

			//		Renderer->AddActor(actor);*/

			//		//ReallyDeletePoint(cutter->GetOutput()->GetPoints(), wrongPointNum);
			//	}

			//	/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			//	sphereSource->SetCenter(pos);
			//	sphereSource->SetRadius(0.1);*/

			//	/*vtkVector3d test = (vtkVector3d)pos;

			//	test = test + insersionVec;

			//	pos[0] = test[0];
			//	pos[1] = test[1];
			//	pos[2] = test[2];*/

			//	vtkSmartPointer<vtkLineSource> lineSource =
			//		vtkSmartPointer<vtkLineSource>::New();
			//	lineSource->SetPoint1(point0);
			//	lineSource->SetPoint2(point1);
			//	lineSource->Update();

			//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			//	mapper->SetInputConnection(lineSource->GetOutputPort());

			//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			//	actor->SetMapper(mapper);
			//	actor->GetProperty()->SetColor(colorRGB(0, 0, 255));
			//	actor->GetProperty()->SetOpacity(0.7);

			//	Renderer->AddActor(actor);

			//	CurrentPointCnt++;

			//	if (CurrentPointCnt > cutterNumberOfPoints)
			//	{
			//		//test++;

			//		cout << "1 cnt : ";

			//		mpPointLocator->SetDataSet(cutter->GetOutput());
			//		mpPointLocator->BuildLocator();

			//		/*if(test==2)
			//		break;*/
			//		if (isFalt)
			//		{
			//			CurrentPointCnt = 0;
			//			isFalt = false;
			//			continue;
			//		}
			//		else
			//			break;
			//	}
			//}

			cout << "delay cnt : " << test;
			//cout << "cutter cnt : " << cutter->GetOutput()->GetPoints()->GetNumberOfPoints();

			//while (1)
			//{
			//	vtkVector3d cutterPoint;
			//	vtkVector3d nextCutterPoint;
			//	vtkVector3d tempPoint1, tempPoint2;

			//	vtkVector3d reverseVec;

			//	reverseVec[0] = 0;
			//	reverseVec[1] = 0;
			//	reverseVec[2] = 2;

			//	cutter->GetOutput()->GetPoints()->GetPoint(test - 1, cutterPoint.GetData());

			//	tempPoint1 = cutterPoint - reverseVec;
			//	tempPoint2 = cutterPoint + reverseVec;

			//	double point0[3];
			//	double point1[3];

			//	point0[0] = tempPoint1[0];
			//	point0[1] = tempPoint1[1];
			//	point0[2] = tempPoint1[2];

			//	point1[0] = tempPoint2[0];
			//	point1[1] = tempPoint2[1];
			//	point1[2] = tempPoint2[2];

			//	mpPointLocator->IntersectWithLine(point0, point1, tolerance, t, pos1, pcoords, subId);

			//	/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			//	sphereSource->SetCenter(pos);
			//	sphereSource->SetRadius(0.1);

			//	vtkSmartPointer<vtkLineSource> lineSource =
			//		vtkSmartPointer<vtkLineSource>::New();
			//	lineSource->SetPoint1(point0);
			//	lineSource->SetPoint2(point1);
			//	lineSource->Update();

			//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			//	mapper->SetInputConnection(sphereSource->GetOutputPort());
			//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			//	actor->SetMapper(mapper);
			//	actor->GetProperty()->SetColor(colorRGB(0, 255, 255));
			//	actor->GetProperty()->SetOpacity(0.7);

			//	Renderer->AddActor(actor);*/

			//	//New UnderLine, 마지막 
			//	cutter->GetOutput()->GetPoints()->SetPoint(test - 1, pos2);

			//	test++;

			//	if (test > cutter->GetOutput()->GetPoints()->GetNumberOfPoints())
			//	{
			//		mpPointLocator->SetDataSet(cutter->GetOutput());
			//		mpPointLocator->BuildLocator();
			//		break;
			//	}
			//}

			//### OutLine Draw###

			//while(1)
			//{
			//	vtkSmartPointer<vtkIdList> CellPointIDList = vtkSmartPointer<vtkIdList>::New();
			//	for (int i = 0; i < 2; i++)
			//	{
			//		cutter->GetOutput()->GetCellPoints(cellIdList->GetId(i), CellPointIDList);

			//		if (CellPointIDList->GetId(0) != NowPointID && CellPointIDList->GetId(0) != PrePointID)
			//		{
			//			PrePointID = NowPointID;
			//			NowPointID = CellPointIDList->GetId(0);
			//			break;
			//		}
			//		else if (CellPointIDList->GetId(1) != NowPointID && CellPointIDList->GetId(1) != PrePointID)
			//		{
			//			PrePointID = NowPointID;
			//			NowPointID = CellPointIDList->GetId(1);
			//			break;
			//		}
			//	}
			//	connectedVertices->InsertNextId(NowPointID);
			//	cutter->GetOutput()->GetPointCells(NowPointID, cellIdList);

			//	ReMakePoints->InsertNextPoint(cutter->GetOutput()->GetPoints()->GetPoint(NowPointID));

			//	vtkSmartPointer<vtkLine> TmpLine = vtkSmartPointer<vtkLine>::New();
			//	TmpLine->GetPointIds()->SetId(0, numCell++);
			//	TmpLine->GetPointIds()->SetId(1, numCell);
			//	ReMakeCells->InsertNextCell(TmpLine);

			//	if (cellIdList->GetNumberOfIds() == 1)
			//		break;
			//}

			//vtkVector3d PreNorVector;
			//std::cerr << "Point Sort Done"<< std::endl;

			///*vtkSmartPointer<vtkPolyData> linesPolyData =
			//	vtkSmartPointer<vtkPolyData>::New();

			//linesPolyData->SetPoints(ReMakePoints);
			//linesPolyData->SetLines(ReMakeCells);


			//{
			//	vtkSmartPointer<vtkPolyDataMapper> Mapper =
			//		vtkSmartPointer<vtkPolyDataMapper>::New();
			//	Mapper->SetInputData(linesPolyData);
			//	Mapper->ScalarVisibilityOff();
			//	vtkSmartPointer<vtkActor> Actor =
			//		vtkSmartPointer<vtkActor>::New();
			//	Actor->SetMapper(Mapper);
			//	Actor->GetProperty()->SetColor(0, 0, 0);
			//	Renderer->AddActor(Actor);
			//}*/

			////vtkVector3d PreNormalVec(NULL, NULL, NULL);
			//for (vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
			//{
			//	vtkVector3d pointA, PointB;
			//	vtkVector3d normalVec, CrossVec, LineToLineVec;

			//	plane->GetNormal(CrossVec.GetData());

			//	cutter->GetOutput()->GetPoints()->GetPoint(connectedVertices->GetId(i), pointA.GetData());

			//	if (i == connectedVertices->GetNumberOfIds() - 1)
			//	{
			//		cutter->GetOutput()->GetPoints()->GetPoint(connectedVertices->GetId(i - 1), PointB.GetData());
			//		vtkMath::Subtract(pointA.GetData(), PointB.GetData(), LineToLineVec.GetData());
			//		break;
			//	}
			//	else
			//	{
			//		cutter->GetOutput()->GetPoints()->GetPoint(connectedVertices->GetId(i + 1), PointB.GetData());
			//		vtkMath::Subtract(PointB.GetData(), pointA.GetData(), LineToLineVec.GetData());
			//	}

			//	vtkMath::Normalize(LineToLineVec.GetData());
			//	//vtkMath::Normalize(CrossVec.GetData());
			//	vtkMath::Cross(LineToLineVec.GetData(), CrossVec.GetData(), normalVec.GetData());

			//	vtkMath::Normalize(normalVec.GetData());

			//	pointA = pointA + normalVec;

			//	if (i != 0)
			//	{
			//		std::cerr << "@ 법선백터와 크로스백터의 각도는? :" << CrossVec.Dot(normalVec) << std::endl;
			//	}
			//	PreNorVector = normalVec;

			//	cutter->GetOutput()->GetPoints()->SetPoint(connectedVertices->GetId(i), pointA.GetData());
			//	cutter->Update();
			//}

			//if (CutterTest != NULL)
			//{
			//	Renderer->RemoveActor(CutterTest);
			//}

			/*CutterTest = vtkSmartPointer<vtkActor>::New();
			vtkSmartPointer<vtkPolyDataMapper> cutterMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			cutterMapper->SetInputConnection(cutter->GetOutputPort());
			cutterMapper->ScalarVisibilityOff();

			CutterTest->GetProperty()->SetColor(1.0, 0, 1.);
			CutterTest->GetProperty()->SetLineWidth(1);
			CutterTest->SetMapper(cutterMapper);
			Renderer->AddActor(CutterTest);
			Renderer->GetRenderWindow()->Render();*/
		}
		else if (key.compare("5") == 0)
		{
			vtkPolyData* sphereData = updateData;

			// Create a data array to hold the weighting coefficients
			vtkSmartPointer<vtkFloatArray> tfarray =
				vtkSmartPointer<vtkFloatArray>::New();
			vtkIdType npoints = sphereData->GetNumberOfPoints();
			tfarray->SetNumberOfComponents(2);
			tfarray->SetNumberOfTuples(npoints);

			// Parameterize the sphere along the z axis, and fill the weights
			// with (1.0-a, a) to linearly interpolate across the shape
			for (int i = 0; i < npoints; i++)
			{
				double pt[3];
				sphereData->GetPoint(i, pt);
				//    double x = pt[0];
				//    double y = pt[1];
				double z = pt[2];

				double zn = z + 1;
				double zn1 = 1.0 - zn;
				if (zn > 1.0)
				{
					zn = 1.0;
				}
				if (zn1 < 0.0)
				{
					zn1 = 0.0;
				}

				tfarray->SetComponent(i, 0, zn1);
				tfarray->SetComponent(i, 1, zn);
			}

			// Create field data to hold the array, and bind it to the sphere
			//  vtkSmartPointer<vtkFieldData> fd =
			//    vtkSmartPointer<vtkFieldData>::New();
			tfarray->SetName("weights");
			sphereData->GetPointData()->AddArray(tfarray);

			// Use an ordinary transform to stretch the shape
			vtkSmartPointer<vtkTransform> stretch =
				vtkSmartPointer<vtkTransform>::New();
			stretch->Scale(1, 1, 1);

			vtkSmartPointer<vtkTransformFilter> stretchFilter =
				vtkSmartPointer<vtkTransformFilter>::New();
#if VTK_MAJOR_VERSION <= 5
			stretchFilter->SetInput(sphereData);
#else
			stretchFilter->SetInputData(sphereData);
#endif
			stretchFilter->SetTransform(stretch);

			// Now, for the weighted transform stuff
			vtkSmartPointer<vtkWeightedTransformFilter> weightedTrans =
				vtkSmartPointer<vtkWeightedTransformFilter>::New();

			// Create two transforms to interpolate between
			vtkSmartPointer<vtkTransform> identity =
				vtkSmartPointer<vtkTransform>::New();
			identity->Identity();

			vtkSmartPointer<vtkTransform> rotated =
				vtkSmartPointer<vtkTransform>::New();
			double rotatedAngle = 0;
			rotated->RotateX(rotatedAngle);

			weightedTrans->SetNumberOfTransforms(2);
			weightedTrans->SetTransform(identity, 0);
			weightedTrans->SetTransform(rotated, 1);
			// which data array should the filter use ?
			weightedTrans->SetWeightArray("weights");

			weightedTrans->SetInputConnection(stretchFilter->GetOutputPort());

			vtkSmartPointer<vtkPolyDataMapper> weightedTransMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			weightedTransMapper->SetInputConnection(weightedTrans->GetOutputPort());
			//weightedTransMapper->SetInputConnection(stretchFilter->GetOutputPort());

			vtkSmartPointer<vtkActor> weightedTransActor =
				vtkSmartPointer<vtkActor>::New();
			weightedTransActor->SetMapper(weightedTransMapper);
			weightedTransActor->GetProperty()->SetDiffuseColor(0.8, 0.8, 0.1);
			weightedTransActor->GetProperty()->SetRepresentationToSurface();

			Renderer->AddActor(weightedTransActor);
			Renderer->GetRenderWindow()->Render();
		}
		else if (key.compare("6") == 0)
		{
			double minBound[3];
			minBound[0] = updateData->GetBounds()[0];
			minBound[1] = updateData->GetBounds()[2];
			minBound[2] = updateData->GetBounds()[4];

			double maxBound[3];
			maxBound[0] = updateData->GetBounds()[1];
			maxBound[1] = updateData->GetBounds()[3];
			maxBound[2] = updateData->GetBounds()[5];

			double center[3];
			center[0] = updateData->GetCenter()[0];
			center[1] = updateData->GetCenter()[1];
			center[2] = updateData->GetCenter()[2];

			//Create an arrow.
			vtkSmartPointer<vtkArrowSource> arrowSource =
				vtkSmartPointer<vtkArrowSource>::New();

			// Generate a random start and end point
			double startPoint[3], endPoint[3];
#ifndef main
			vtkMath::RandomSeed(time(NULL));
#else
			vtkMath::RandomSeed(8775070);
#endif
			/*startPoint[0] = vtkMath::Random(-10, 10);
			startPoint[1] = vtkMath::Random(-10, 10);
			startPoint[2] = vtkMath::Random(-10, 10);
			endPoint[0] = vtkMath::Random(-10, 10);
			endPoint[1] = vtkMath::Random(-10, 10);
			endPoint[2] = vtkMath::Random(-10, 10);*/

			startPoint[0] = center[0];
			startPoint[1] = center[1] / 4;
			startPoint[2] = maxBound[2];
			endPoint[0] = center[0];
			endPoint[1] = center[1];
			endPoint[2] = center[2];

			// Compute a basis
			double normalizedX[3];
			double normalizedY[3];
			double normalizedZ[3];

			// The X axis is a vector from start to end
			vtkMath::Subtract(endPoint, startPoint, normalizedX);
			double length = vtkMath::Norm(normalizedX);
			vtkMath::Normalize(normalizedX);

			// The Z axis is an arbitrary vector cross X
			double arbitrary[3];
			arbitrary[0] = vtkMath::Random(-10, 10);
			arbitrary[1] = vtkMath::Random(-10, 10);
			arbitrary[2] = vtkMath::Random(-10, 10);
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
			transformPD->SetInputConnection(arrowSource->GetOutputPort());

			//Create a mapper and actor for the arrow
			vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
#ifdef USER_MATRIX
			mapper->SetInputConnection(arrowSource->GetOutputPort());
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
			Renderer->AddActor(actor);
			Renderer->AddActor(sphereStart);
			Renderer->AddActor(sphereEnd);
			Renderer->SetBackground(.1, .2, .3); // Background color dark blue

			//Render and interact
			/*renderWindow->Render();
			renderWindowInteractor->Start();*/
		}
		else if (key.compare("7") == 0)
		{
			vtkSmartPointer<vtkPolyDataMapper> mapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
			mapper->SetInput(sphere);
#else
			mapper->SetInputData(updateData);
#endif
			vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);

			// Create the outline
			vtkSmartPointer<vtkOutlineFilter> outline =
				vtkSmartPointer<vtkOutlineFilter>::New();
#if VTK_MAJOR_VERSION <= 5
			outline->SetInput(sphere);
#else
			outline->SetInputData(updateData);
#endif
			vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
				vtkSmartPointer<vtkPolyDataMapper>::New();
			outlineMapper->SetInputConnection(outline->GetOutputPort());
			vtkSmartPointer<vtkActor> outlineActor =
				vtkSmartPointer<vtkActor>::New();
			outlineActor->SetMapper(outlineMapper);
			outlineActor->GetProperty()->SetColor(0, 0, 0);

			// Setup the window
			vtkSmartPointer<vtkRenderer> renderer =
				vtkSmartPointer<vtkRenderer>::New();
			vtkSmartPointer<vtkRenderWindow> renderWindow =
				vtkSmartPointer<vtkRenderWindow>::New();
			renderWindow->AddRenderer(renderer);
			vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
				vtkSmartPointer<vtkRenderWindowInteractor>::New();
			renderWindowInteractor->SetRenderWindow(renderWindow);

			// Add the actors to the scene
			Renderer->AddActor(actor);
			Renderer->AddActor(outlineActor);
			Renderer->SetBackground(1, 1, 1); // Background color white
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

		//vtkSmartPointer<vtkPolyDataMapper> mapper =
		//	vtkSmartPointer<vtkPolyDataMapper>::New();
		////inputMapper->SetInputData(input);
		//mapper->SetInputData(updateData);
		//mapper->ScalarVisibilityOff();
		////vtkSmartPointer<vtkActor> Actor =			vtkSmartPointer<vtkActor>::New();
		//Actor->SetMapper(mapper);
		//Actor->GetProperty()->SetSpecular(0.5);
		//Actor->GetProperty()->SetDiffuse(.2);
		//Actor->GetProperty()->SetAmbient(1.0);
		//Actor->GetProperty()->SetAmbientColor(.4, .4, .0);
		//Actor->GetProperty()->SetColor(.5, .5, .2); 
		//vtkSmartPointer<vtkProperty> backFaces =
		//	vtkSmartPointer<vtkProperty>::New();
		//backFaces->SetSpecular(0.5);
		//backFaces->SetDiffuse(.2);
		//backFaces->SetAmbient(1.0);
		//backFaces->SetAmbientColor(.4, .4, .0);
		//backFaces->SetColor(.2, .2, .6);
		//Actor->SetBackfaceProperty(backFaces);
		////Renderer->AddActor(Actor);

		//if (value)
		//	Renderer->AddActor(inputActor);
		//else
		//	Renderer->RemoveActor(inputActor);
		//Renderer->GetRenderWindow()->Render();
		//Renderer->Render();
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

	int CompareDoubleAbsoulte(double *x, double *y)
	{
		double absTolerance = 0.01;

		double diffX = x[0] - y[0];
		double diffY = x[1] - y[1];
		double diffZ = x[2] - y[2];

		if (diffX < -absTolerance || diffX > absTolerance || diffY < -absTolerance || diffY > absTolerance || diffZ < -absTolerance || diffZ > absTolerance)
			return 0;
		else
			return 1;
	}

	/*void ReallyDeletePoint(vtkSmartPointer<vtkPoints> points, vtkIdType id)
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
	}*/

	vtkSmartPointer<vtkActor> outActor = NULL;
	vtkSmartPointer<vtkActor> CutterTest = NULL;;
	vtkActor* Actor;
	vtkActor* inputActor;
	vtkPolyData* updateData;
	vtkPolyData* inputData;
	vtkPolyData* BooleanTargetData;
	vtkPolyData* SleeveData;
	vtkRenderer* Renderer;
};
vtkStandardNewMacro(KeyPressInteractorStyle);


vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);

	/*
	cout << "Vertex 0 is used in cells ";
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
	cout << cellIdList->GetId(i) << ", ";
	}
	cout << endl;
	*/

	for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if (pointIdList->GetId(0) != id)
		{
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			//cout << "Connected to " << pointIdList->GetId(1) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}

	return connectedVertices;
}

bool GetPointNormals(vtkPolyData* polydata)
{
	std::cout << "In GetPointNormals: " << polydata->GetNumberOfPoints() << std::endl;
	std::cout << "Looking for point normals..." << std::endl;

	// Count points
	vtkIdType numPoints = polydata->GetNumberOfPoints();
	std::cout << "There are " << numPoints << " points." << std::endl;

	// Count triangles
	vtkIdType numPolys = polydata->GetNumberOfPolys();
	std::cout << "There are " << numPolys << " polys." << std::endl;

	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkDoubleArray* normalDataDouble =
		vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

	if (normalDataDouble)
	{
		int nc = normalDataDouble->GetNumberOfTuples();
		std::cout << "There are " << nc
			<< " components in normalDataDouble" << std::endl;
		return true;
	}

	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkFloatArray* normalDataFloat =
		vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

	if (normalDataFloat)
	{
		int nc = normalDataFloat->GetNumberOfTuples();
		std::cout << "There are " << nc
			<< " components in normalDataFloat" << std::endl;
		return true;
	}

	////////////////////////////////////////////////////////////////
	// Point normals
	vtkDoubleArray* normalsDouble =
		vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetNormals());

	if (normalsDouble)
	{
		std::cout << "There are " << normalsDouble->GetNumberOfComponents()
			<< " components in normalsDouble" << std::endl;
		return true;
	}

	////////////////////////////////////////////////////////////////
	// Point normals
	vtkFloatArray* normalsFloat =
		vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetNormals());

	if (normalsFloat)
	{
		std::cout << "There are " << normalsFloat->GetNumberOfComponents()
			<< " components in normalsFloat" << std::endl;
		return true;
	}

	/////////////////////////////////////////////////////////////////////
	// Generic type point normals
	vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals(); //works
	if (normalsGeneric)
	{
		std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
			<< " normals in normalsGeneric" << std::endl;

		double testDouble[3];
		normalsGeneric->GetTuple(0, testDouble);

		std::cout << "Double: " << testDouble[0] << " "
			<< testDouble[1] << " " << testDouble[2] << std::endl;

		// Can't do this:
		/*
		float testFloat[3];
		normalsGeneric->GetTuple(0, testFloat);

		std::cout << "Float: " << testFloat[0] << " "
		<< testFloat[1] << " " << testFloat[2] << std::endl;
		*/
		return true;
	}


	// If the function has not yet quit, there were none of these types of normals
	std::cout << "Normals not found!" << std::endl;
	return false;

}

int main(int, char *[])
{
	//vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszPathName2 = { "D:\\stl\\sample.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();

	const char* lpszPathName3 = { "D:\\stl\\sample.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkSTLReader> reader2 =
		vtkSmartPointer<vtkSTLReader>::New();
	reader2->SetFileName(lpszPathName3);
	reader2->Update();

	vtkSmartPointer<vtkPolyData> Lineinput = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszLinePathName = { "D:\\TestPoints.vtp" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkXMLPolyDataReader> Linereader =
		vtkSmartPointer<vtkXMLPolyDataReader>::New();
	Linereader->SetFileName(lpszLinePathName);
	Linereader->Update();

	input->DeepCopy(reader->GetOutput());
	input2->DeepCopy(reader2->GetOutput());
	//Lineinput->DeepCopy(Linereader->GetOutput()->GetPoints());
	vtkSmartPointer<vtkPoints> points = Linereader->GetOutput()->GetPoints();
	//points->Modified();

	vtkSmartPointer<vtkParametricSpline>mpSelectSpline = vtkSmartPointer<vtkParametricSpline>::New();
	mpSelectSpline->SetPoints(points);
	mpSelectSpline->Modified();
	vtkSmartPointer<vtkParametricFunctionSource> FunctionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
	FunctionSource->SetParametricFunction(mpSelectSpline);
	FunctionSource->SetUResolution(500);
	FunctionSource->Update();

	std::cout << "Before decimation" << std::endl << "------------" << std::endl;
	std::cout << "There are " << input->GetNumberOfPoints() << " points." << std::endl;
	std::cout << "There are " << input->GetNumberOfPolys() << " polygons." << std::endl;
	std::cout << "There are " << input->GetNumberOfLines() << " polygons." << std::endl;

	vtkSmartPointer<vtkPolyDataMapper> lineMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	lineMapper->SetInputData(FunctionSource->GetOutput());
	//lineMapper->SetlineConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> lineActor =
		vtkSmartPointer<vtkActor>::New();
	lineActor->SetMapper(lineMapper);
	lineActor->GetProperty()->SetColor(1., 0., 0.);

	vtkSmartPointer<vtkPolyDataMapper> inputMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	inputMapper->SetInputData(input);


	vtkSmartPointer<vtkTransform> translation =
		vtkSmartPointer<vtkTransform>::New();
	translation->Translate(1.0, 1.0, 1.0);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetInputData(input2);
	transformFilter->SetTransform(translation);
	transformFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> inputMapper2 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	//inputMapper2->SetInputData(transformFilter->GetOutput());
	inputMapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);


	vtkSmartPointer<vtkActor> inputActor2 =
		vtkSmartPointer<vtkActor>::New();
	inputActor2->SetMapper(inputMapper2);
	inputActor2->GetProperty()->SetColor(1., 0., 0.);

	vtkSmartPointer<vtkPolyDataMapper> decimatedMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	decimatedMapper->SetInputConnection(reader->GetOutputPort());
	vtkSmartPointer<vtkActor> decimatedActor =
		vtkSmartPointer<vtkActor>::New();
	decimatedActor->SetMapper(decimatedMapper);
	decimatedActor->GetProperty()->SetColor(.2, .2, .2);

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
	rightRenderer->SetBackground(1, 1, 1);
	rightRenderer->SetActiveCamera(camera);

	// Add the sphere to the left and the cube to the right
	//leftRenderer->AddActor(lineActor);
	//leftRenderer->AddActor(inputActor);
	//leftRenderer->AddActor(inputActor2);
	//rightRenderer->AddActor(lineActor);
	//rightRenderer->AddActor(decimatedActor);
	rightRenderer->AddActor(inputActor);
	//rightRenderer->AddActor(surfaceActor);

	//leftRenderer->AddActor(contoursActor);
	//leftRenderer->AddActor(surfaceActor);

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
	std::cout << "There are " << update->GetNumberOfPoints() << " points." << std::endl;
	style->updateData = update;
	style->Renderer = rightRenderer;
	interactor->SetInteractorStyle(style);
	style->SetCurrentRenderer(rightRenderer);

	//// Create Axis Widget
	vtkSmartPointer<vtkAxesActor>mpMainAxes = vtkSmartPointer<vtkAxesActor>::New();
	vtkSmartPointer<vtkOrientationMarkerWidget> mpMainAxisWidget =
		vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	//mpMainAxisWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
	mpMainAxisWidget->SetOrientationMarker(mpMainAxes);
	mpMainAxisWidget->SetInteractor(interactor);
	mpMainAxisWidget->SetViewport(0.9, 0.0, 1.0, 0.1);
	mpMainAxisWidget->SetEnabled(1);
	mpMainAxisWidget->InteractiveOff();

	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}


