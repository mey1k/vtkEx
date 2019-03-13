else if (key.compare("4") == 0)
{
	vtkSmartPointer<vtkPlane> plane =
		vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(updateData->GetCenter());
	//plane->SetOrigin(10,10,5);
	plane->SetNormal(-0.5, 1, 0);

	//plane->SetNormal(1, 1, 1);

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

	// Create cutter
	vtkSmartPointer<vtkCutter> cutter =
		vtkSmartPointer<vtkCutter>::New();
	cutter->SetCutFunction(plane);
	cutter->SetInputData(updateData);
	cutter->GenerateValues(200, -distanceMin, distanceMax);
	cutter->Update();

	vtkSmartPointer<vtkCutter> copyCutter = vtkSmartPointer<vtkCutter>::New();
	copyCutter->SetCutFunction(plane);
	copyCutter->SetInputData(updateData);
	cutter->GenerateValues(200, -distanceMin, distanceMax);
	copyCutter->Update();

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
		cellIdList = vtkSmartPointer<vtkIdList>::New();
		cutter->GetOutput()->GetPointCells(NowPointID, cellIdList);

		if (cellIdList->GetNumberOfIds() == 1)
		{
			//std::cerr << "@NowPointID :" << NowPointID <<" | Start or end"<< std::endl;
			std::cerr << "@cellID:" << cellIdList->GetId(0) << " | Start or end CELL" << std::endl;
			connectedVertices->InsertNextId(NowPointID);
			PrePointID = NowPointID;
			break;
		}
	}

	vtkSmartPointer<vtkPoints> ReMakePoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> ReMakeCells = vtkSmartPointer<vtkCellArray>::New();
	int numCell = 0;
	ReMakePoints->InsertNextPoint(cutter->GetOutput()->GetPoints()->GetPoint(NowPointID));

	// Create the locator
	/*vtkSmartPointer<vtkCellLocator> mpPointLocator =
	vtkSmartPointer<vtkCellLocator>::New();
	mpPointLocator->SetDataSet(cutter->GetOutput());
	mpPointLocator->BuildLocator();*/

	vtkSmartPointer<vtkCellLocator> mpPointLocator =
		vtkSmartPointer<vtkCellLocator>::New();
	mpPointLocator->SetDataSet(cutter->GetOutput());
	mpPointLocator->BuildLocator();

	int cutterNumberOfPoints = cutter->GetOutput()->GetPoints()->GetNumberOfPoints();

	double tolerance = 0.001;
	double pos1[] = { 0.0, 0.0, 0.0 };
	double pos2[] = { 0.0, 0.0, 0.0 };
	double pcoords[] = { 0.0, 0.0, 0.0 };
	int subId = 0;
	double t = 0;
	int CurrentPointCnt = 1;

	vtkVector3d cutterDirection, startPoint, endPoint;
	//start
	cutter->GetOutput()->GetPoints()->GetPoint(0, startPoint.GetData());
	//end
	cutter->GetOutput()->GetPoints()->GetPoint(cutterNumberOfPoints - 1, endPoint.GetData());
	//direction
	vtkMath::Subtract(endPoint.GetData(), startPoint.GetData(), cutterDirection.GetData());
	vtkMath::Normalize(cutterDirection.GetData());

	vtkVector3d insersionVec;
	insersionVec[0] = 0;
	insersionVec[1] = 0;
	insersionVec[2] = 5;

	vtkVector3d cutterPoint;
	vtkVector3d topPoint, bottomPoint;

	bool isFalt = false;

	double dValue = 1.05;
	double value = floor(100.*(dValue + 0.005)) / 100.;
	int test = 0;


	double bounds[6];

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




	while (1)
	{
		copyCutter->GetOutput()->GetPoints()->GetPoint(CurrentPointCnt - 1, cutterPoint.GetData());

		topPoint = cutterPoint + insersionVec;
		bottomPoint = cutterPoint - insersionVec;

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

		//new Line
		if (!CompareDoubleAbsoulte(pos1, pos2))
		{
			isFalt = true;

			int wrongPointNum = cutter->GetOutput()->FindPoint(pos2[0], pos2[1], pos2[2]);


			vtkVector3d wrongPointVec = (vtkVector3d)pos2;
			vtkVector3d modifiedVec;

			if (wrongPointNum < cutterNumberOfPoints / 2)
			{
				modifiedVec = (wrongPointVec + cutterDirection) - cutterDirection * value;
			}
			else
			{
				modifiedVec = (wrongPointVec - cutterDirection) + cutterDirection * value;
			}

			pos2[0] = modifiedVec[0];
			pos2[1] = modifiedVec[1];
			pos2[2] = modifiedVec[2];

			cutter->GetOutput()->GetPoints()->SetPoint(wrongPointNum, pos2);


			/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->SetCenter(pos2);
			sphereSource->SetRadius(0.1);

			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(sphereSource->GetOutputPort());

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(colorRGB(0, 0, 255));
			actor->GetProperty()->SetOpacity(0.7);

			Renderer->AddActor(actor);*/

			//ReallyDeletePoint(cutter->GetOutput()->GetPoints(), wrongPointNum);
		}

		/*vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter(pos);
		sphereSource->SetRadius(0.1);*/

		/*vtkVector3d test = (vtkVector3d)pos;

		test = test + insersionVec;

		pos[0] = test[0];
		pos[1] = test[1];
		pos[2] = test[2];*/

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

		CurrentPointCnt++;

		if (CurrentPointCnt > cutterNumberOfPoints)
		{
			//test++;

			cout << "1 cnt : ";

			mpPointLocator->SetDataSet(cutter->GetOutput());
			mpPointLocator->BuildLocator();

			/*if(test==2)
			break;*/
			if (isFalt)
			{
				CurrentPointCnt = 0;
				isFalt = false;
				continue;
			}
			else
				break;
		}
	}

	cout << "delay cnt : " << test;
	cout << "cutter cnt : " << cutter->GetOutput()->GetPoints()->GetNumberOfPoints();

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

	CutterTest = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkPolyDataMapper> cutterMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cutterMapper->SetInputConnection(cutter->GetOutputPort());
	cutterMapper->ScalarVisibilityOff();

	CutterTest->GetProperty()->SetColor(1.0, 0, 1.);
	CutterTest->GetProperty()->SetLineWidth(1);
	CutterTest->SetMapper(cutterMapper);
	Renderer->AddActor(CutterTest);
	Renderer->GetRenderWindow()->Render();
}