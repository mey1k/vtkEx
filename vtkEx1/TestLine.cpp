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
//#include <vtkCollisionDetectionFilter>

vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();

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

			// Create cutter
			vtkSmartPointer<vtkCutter> cutter =
				vtkSmartPointer<vtkCutter>::New();
			cutter->SetCutFunction(plane);
			//cutter->SetSortByToSortByValue();
			//cutter->SetSortByToSortByCell();
			cutter->SetInputData(updateData);
			//cutter->GenerateValues(10, -distanceMin, distanceMax);
			cutter->Update();


			//vtkPolyData* cutterCopy = vtkPolyData::New();
			//cutterCopy->DeepCopy(cutter->GetOutput());


			vtkSmartPointer<vtkPolyDataNormals> normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
			normalFilter->SetInputData(cutter->GetOutput());
			normalFilter->ComputePointNormalsOn();
			normalFilter->Update();
			normalFilter->GetOutput()->GetPoint((vtkIdType)1);
			
			vtkIdType NowPointID, PrePointID;
			vtkSmartPointer<vtkIdList> cellIdList;
			vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

			for (NowPointID = 0; NowPointID < cutter->GetOutput()->GetPoints()->GetNumberOfPoints(); NowPointID++)
			{
				cellIdList =	vtkSmartPointer<vtkIdList>::New();
				cutter->GetOutput()->GetPointCells(NowPointID, cellIdList);
				if (cellIdList->GetNumberOfIds() == 1)
				{
					//std::cerr << "@NowPointID :" << NowPointID <<" | Start or end"<< std::endl;
					std::cerr << "@cellID:" << cellIdList ->GetId(0)<<" | Start or end CELL"<< std::endl;
					connectedVertices->InsertNextId(NowPointID);
					PrePointID = NowPointID;
					break;
				}
			}

			vtkSmartPointer<vtkPoints> ReMakePoints = vtkSmartPointer<vtkPoints>::New();
			vtkSmartPointer<vtkCellArray> ReMakeCells = vtkSmartPointer<vtkCellArray>::New();
			int numCell = 0;
			ReMakePoints->InsertNextPoint(cutter->GetOutput()->GetPoints()->GetPoint(NowPointID));


			//kmj

			//vtkSmartPointer<vtkPolyData> cutterPolydata =
			//	vtkSmartPointer<vtkPolyData>::New();

			//cutterPolydata->SetPoints(cutter->GetOutput()->GetPoints());


			vtkSmartPointer<vtkPolyLine> cutterLine = vtkSmartPointer<vtkPolyLine>::New();
			vtkSmartPointer<vtkCellArray> cutterCell = vtkSmartPointer<vtkCellArray>::New();
			vtkSmartPointer<vtkPolyData> cutterPolyData = vtkSmartPointer<vtkPolyData>::New();

			//cutterLine->Initialize();
			//cutterLine->GetPointIds()->SetNumberOfIds(cutter->GetOutput()->GetPoints()->GetNumberOfPoints());
			//for (unsigned int i = 0; i < cutter->GetOutput()->GetPoints()->GetNumberOfPoints(); i++)
			//{
			//	cutterLine->GetPointIds()->SetId(i, i);
			//}
			//cutterLine->Modified();

			//cutterCell->Initialize();
			//cutterCell->InsertNextCell(cutterLine);
			//cutterCell->Modified();

			//cutterPolyData->Initialize();
			////cutter->GetOutput()->GetPointData
			//cutterPolyData->SetPoints(cutter->GetOutput()->GetPoints());
			//cutterPolyData->SetLines(cutterCell);
			//cutterPolyData->Modified();

			vtkSmartPointer<vtkCellLocator>	mpPointLocator = vtkSmartPointer<vtkCellLocator>::New();
			mpPointLocator->SetDataSet(cutter->GetOutput());
			mpPointLocator->BuildLocator();

			/*vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputData(cutterPolyData);
			mapper->Modified();
			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(colorRGB(255, 0, 255));
			actor->GetProperty()->SetLineWidth(1);

			Renderer->AddActor(actor);*/

			double tolerance = 0.001;
			double pos[] = { 0.0, 0.0, 0.0 };
			double pcoords[] = { 0.0, 0.0, 0.0 };
			int subId = 0;
			double t = 0;
			int CurrentPointCnt = 1;

			int test = 1;

			while (1)
			{
				vtkVector3d cutterPoint;
				vtkVector3d nextCutterPoint;
				vtkVector3d tempPoint;

				vtkVector3d reverseVec;

				reverseVec[0] = 0;
				reverseVec[1] = 0;
				reverseVec[2] = 2;

				cutter->GetOutput()->GetPoints()->GetPoint(CurrentPointCnt - 1, cutterPoint.GetData());

				tempPoint = cutterPoint + reverseVec;

				double point0[3];
				double point1[3];

				point0[0] = cutterPoint[0];
				point0[1] = cutterPoint[1];
				point0[2] = cutterPoint[2];

				point1[0] = tempPoint[0];
				point1[1] = tempPoint[1];
				point1[2] = tempPoint[2];

				mpPointLocator->IntersectWithLine(point0, point1, tolerance, t, pos, pcoords, subId);

				/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
				sphereSource->SetCenter(pos);
				sphereSource->SetRadius(0.1);*/

				vtkVector3d test = (vtkVector3d)pos;

				test = test + reverseVec;

				pos[0] = test[0];
				pos[1] = test[1];
				pos[2] = test[2];

				/*vtkSmartPointer<vtkLineSource> lineSource =
				vtkSmartPointer<vtkLineSource>::New();
				lineSource->SetPoint1(point0);
				lineSource->SetPoint2(point1);
				lineSource->Update();

				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputConnection(sphereSource->GetOutputPort());
				
				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(colorRGB(0, 0, 255));
				actor->GetProperty()->SetOpacity(0.7);

				Renderer->AddActor(actor);*/

				//New UnderLine, 마지막 
				cutter->GetOutput()->GetPoints()->SetPoint(CurrentPointCnt - 1, pos);

				CurrentPointCnt++;

				if (CurrentPointCnt > cutter->GetOutput()->GetPoints()->GetNumberOfPoints())
				{
					mpPointLocator->SetDataSet(cutter->GetOutput());
					mpPointLocator->BuildLocator();

					break;
				}
			}

			cout << "3ms delay";

			while (1)
			{
				vtkVector3d cutterPoint;
				vtkVector3d nextCutterPoint;
				vtkVector3d tempPoint1, tempPoint2;

				vtkVector3d reverseVec;

				reverseVec[0] = 0;
				reverseVec[1] = 0;
				reverseVec[2] = 2;

				cutter->GetOutput()->GetPoints()->GetPoint(test - 1, cutterPoint.GetData());

				tempPoint1 = cutterPoint - reverseVec;
				tempPoint2 = cutterPoint + reverseVec;

				double point0[3];
				double point1[3];

				point0[0] = tempPoint1[0];
				point0[1] = tempPoint1[1];
				point0[2] = tempPoint1[2];

				point1[0] = tempPoint2[0];
				point1[1] = tempPoint2[1];
				point1[2] = tempPoint2[2];

				mpPointLocator->IntersectWithLine(point0, point1, tolerance, t, pos, pcoords, subId);

				/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
				sphereSource->SetCenter(pos);
				sphereSource->SetRadius(0.1);

				vtkSmartPointer<vtkLineSource> lineSource =
					vtkSmartPointer<vtkLineSource>::New();
				lineSource->SetPoint1(point0);
				lineSource->SetPoint2(point1);
				lineSource->Update();

				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputConnection(sphereSource->GetOutputPort());
				vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
				actor->SetMapper(mapper);
				actor->GetProperty()->SetColor(colorRGB(0, 255, 255));
				actor->GetProperty()->SetOpacity(0.7);

				Renderer->AddActor(actor);*/

				//New UnderLine, 마지막 
				cutter->GetOutput()->GetPoints()->SetPoint(test - 1, pos);

				test++;

				if (test > cutter->GetOutput()->GetPoints()->GetNumberOfPoints())
				{
					mpPointLocator->SetDataSet(cutter->GetOutput());
					mpPointLocator->BuildLocator();
					break;
				}
			}

			cout << "재검출 수정 3ms";

			//### OutLine Draw###
			//while(1)
			//{
			//	bool FindNextCell = false;
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

			//CutterTest = vtkSmartPointer<vtkActor>::New();
			//vtkSmartPointer<vtkPolyDataMapper> cutterMapper =
			//	vtkSmartPointer<vtkPolyDataMapper>::New();
			//cutterMapper->SetInputConnection(cutter->GetOutputPort());
			//cutterMapper->ScalarVisibilityOff();

			//CutterTest->GetProperty()->SetColor(1.0, 0, 1.);
			//CutterTest->GetProperty()->SetLineWidth(1);
			//CutterTest->SetMapper(cutterMapper);
			//Renderer->AddActor(CutterTest);
			//Renderer->GetRenderWindow()->Render();
		}
		else if (key.compare("5") == 0)
		{
		}
		else if (key.compare("6") == 0)
		{
		}
		else if (key.compare("7") == 0)
		{
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
		Actor->GetProperty()->SetSpecular(0.5);
		Actor->GetProperty()->SetDiffuse(.2);
		Actor->GetProperty()->SetAmbient(1.0);
		Actor->GetProperty()->SetAmbientColor(.4, .4, .0);
		Actor->GetProperty()->SetColor(.5, .5, .2); 
		vtkSmartPointer<vtkProperty> backFaces =
			vtkSmartPointer<vtkProperty>::New();
		backFaces->SetSpecular(0.5);
		backFaces->SetDiffuse(.2);
		backFaces->SetAmbient(1.0);
		backFaces->SetAmbientColor(.4, .4, .0);
		backFaces->SetColor(.2, .2, .6);
		Actor->SetBackfaceProperty(backFaces);
		//Renderer->AddActor(Actor);

		if (value)
			Renderer->AddActor(inputActor);
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

int main(int, char *[])
{
	//vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszPathName2 = { "D:\\target1.stl" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(lpszPathName2);
	reader->Update();

	vtkSmartPointer<vtkPolyData> Lineinput = vtkSmartPointer<vtkPolyData>::New();
	const char* lpszLinePathName = { "D:\\TestPoints.vtp" };//{ "C:\\Users\\SeongHun\\OneDrive\\Project\\cowHead.vtp" };
	vtkSmartPointer<vtkXMLPolyDataReader> Linereader =
		vtkSmartPointer<vtkXMLPolyDataReader>::New();
	Linereader->SetFileName(lpszLinePathName);
	Linereader->Update();

	input->DeepCopy(reader->GetOutput());
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
	//inputMapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> inputActor =
		vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);

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
	leftRenderer->AddActor(lineActor);
	leftRenderer->AddActor(inputActor);
	//rightRenderer->AddActor(lineActor);
	rightRenderer->AddActor(decimatedActor);
	rightRenderer->AddActor(inputActor);

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


