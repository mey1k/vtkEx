//#include "stdafx.h"
//#include "vtkAutoInit.h" 
//#include "stdafx.h"
//#include <vtkSmartPointer.h>
//#include <vtkProperty.h>
//#include <vtkSelectPolyData.h>
//#include <vtkSphereSource.h>
//#include <vtkClipPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkLODActor.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkXMLPolyDataReader.h>
//#include <vtkNew.h>
//#include <vtkParametricFunctionSource.h>
//
//VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
//VTK_MODULE_INIT(vtkInteractionStyle);
//
//#define vtkRenderingCore_AUTOINIT 2(vtkRenderingOpenGL2, vtkInteractionStyle)
//#define vtkNew(type,name) \
//  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
//
//int main(int, char *[])
//{
//	
//	const char* lpszPathName1 = { "D:\\testPolyData.vtp" };
//
//	vtkSmartPointer<vtkXMLPolyDataReader> reader1 =
//		vtkSmartPointer<vtkXMLPolyDataReader>::New();
//	reader1->SetFileName(lpszPathName1);
//	reader1->Update();
//
//	const char* lpszPathName2 = { "D:\\testPointData.vtp" };
//
//	vtkSmartPointer<vtkXMLPolyDataReader> reader2 =
//		vtkSmartPointer<vtkXMLPolyDataReader>::New();
//	reader2->SetFileName(lpszPathName2);
//	reader2->Update();
//
//	vtkPolyData* polydata = reader1->GetOutput();
//	vtkPolyData* pointdata = reader2->GetOutput();
//
//	///////////////////////////////////////
//	vtkNew(vtkPolyData, TmpPolydata);
//
//	/* Original Spline*/
//	vtkNew(vtkParametricFunctionSource, FunctionSource);
//	FunctionSource->SetParametricFunction(SelectSpline);
//	FunctionSource->SetUResolution(500);
//	FunctionSource->Update();
//
//	vtkNew(vtkCleanPolyData, cp);
//
//	cp->SetInputData(Polydata);
//	cp->Update();
//
//	vtkPolyDataNormals* norma = vtkPolyDataNormals::New();
//	norma->SetInputData(cp->GetOutput());
//	norma->ComputeCellNormalsOff();
//	norma->ComputePointNormalsOff();
//	norma->Update();
//
//	vtkNew(vtkSelectPolyData, mpSelectedLoopArea);
//	mpSelectedLoopArea->SetInputData(norma->GetOutput());
//	mpSelectedLoopArea->SetLoop(FunctionSource->GetOutput()->GetPoints()); // 루프 설정할 포인터 넘김
//	mpSelectedLoopArea->GenerateSelectionScalarsOn();
//	mpSelectedLoopArea->SetSelectionModeToClosestPointRegion();
//	//mpSelectedLoopArea->SetSelectionModeToSmallestRegion();
//	mpSelectedLoopArea->Update();
//
//	vtkNew(vtkClipPolyData, clip);
//	clip->SetInputConnection(mpSelectedLoopArea->GetOutputPort());
//	clip->Update();
//
//	vtkNew(vtkPolyDataConnectivityFilter, confilter);
//	confilter->SetInputConnection(clip->GetOutputPort());
//	confilter->SetExtractionModeToLargestRegion();
//	confilter->Update();
//
//	vtkNew(vtkPolyDataNormals, clipNormalFilter);
//	clipNormalFilter->SetInputData(confilter->GetOutput());
//	clipNormalFilter->Update();
//
//	TmpPolydata->DeepCopy(clipNormalFilter->GetOutput());
//
//	vtkDataArray* ClipDataNormalArray = TmpPolydata->GetPointData()->GetNormals();
//
//	std::cerr << "number of UpdataData Points :" << TmpPolydata->GetPoints()->GetNumberOfPoints() << std::endl;
//	std::cerr << "number of UpdataData Normal :" << ClipDataNormalArray->GetNumberOfTuples() << std::endl;
//
//	for (int i = 0; i < TmpPolydata->GetPoints()->GetNumberOfPoints(); i++)
//	{
//		//vtkP(vtkSphereSource, testSource);
//		vtkVector3d Point(TmpPolydata->GetPoints()->GetPoint(i));
//		vtkVector3d NormalV(ClipDataNormalArray->GetTuple(i));
//
//		NormalV = NormalV * .4;
//		Point = Point + NormalV;
//
//		TmpPolydata->GetPoints()->SetPoint(i, Point.GetData());
//	}
//
//	/*vtkSmartPointer<vtkSTLWriter> stlWriter =
//	vtkSmartPointer<vtkSTLWriter>::New();
//	const char* lpszPathName = { "D:\\Test.stl" };
//	stlWriter->SetFileName(lpszPathName);
//	stlWriter->SetInputConnection(clipNormalFilter->GetOutputPort());
//	stlWriter->Write();*/
//
//	//vtkBooleanOperationPolyDataFilter* Sleeve_BooleanOperation =
//	//	vtkBooleanOperationPolyDataFilter::New();
//	////	booleanOperation->SetOperationToUnion();
//	////	booleanOperation->SetOperationToIntersection();
//	//Sleeve_BooleanOperation->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
//	//Sleeve_BooleanOperation->SetInputData(0, mpSleevePolydata);
//	//Sleeve_BooleanOperation->SetInputData(1, clipNormalFilter->GetOutput());
//	//Sleeve_BooleanOperation->Update();
//
//	//CheckModel(TmpPolydata);
//
//	///* Step 1. */
//
//	vtkImplicitModeller* implicitModeller =
//		vtkImplicitModeller::New();
//	implicitModeller->SetInputData(TmpPolydata);
//	implicitModeller->AdjustBoundsOn();
//	implicitModeller->SetAdjustDistance(.1); // Adjust by 10%
//	implicitModeller->SetMaximumDistance(.05);
//	implicitModeller->Update();
//
//	double bounds[6];
//	implicitModeller->GetOutput()->GetBounds(bounds);
//	double xrange = bounds[1] - bounds[0];
//
//	vtkContourFilter* contourFilter = vtkContourFilter::New();
//	contourFilter->SetInputConnection(implicitModeller->GetOutputPort());
//	contourFilter->GenerateTrianglesOn();
//	contourFilter->SetValue(0, xrange / 30.0); // 30% of xrange
//	contourFilter->Update();
//
//	vtkReverseSense* reverseIN = vtkReverseSense::New();;
//	reverseIN->SetInputConnection(contourFilter->GetOutputPort());
//	reverseIN->ReverseNormalsOn();
//	reverseIN->Update();
//
//	vtkPolyDataNormals* NormalFilter = vtkPolyDataNormals::New();
//	NormalFilter->SetInputData(reverseIN->GetOutput());
//	NormalFilter->Update();
//	//CheckModel(NormalFilter->GetOutput());
//	///
//
//
//	///* Step 2. */
//	vtkBox* BodyClipBound = vtkBox::New();
//
//	double BodyBounds[6];
//	NormalFilter->GetOutput()->GetBounds(BodyBounds);
//
//	//std::cerr << "Print 6 Value before : "
//	//	<< BodyBounds[0] << ","
//	//	<< BodyBounds[1] << "," 
//	//	<< BodyBounds[2] << "," 
//	//	<< BodyBounds[3] << "," 
//	//	<< BodyBounds[4] << "," 
//	//	<< BodyBounds[5] <<  std::endl;
//
//	//for (int i = 0; i < 6; i++)
//	//	BodyBounds[i] *= 1.5;
//
//	//std::cerr << "Print 6 Value after : "
//	//	<< BodyBounds[0] << ","
//	//	<< BodyBounds[1] << ","
//	//	<< BodyBounds[2] << ","
//	//	<< BodyBounds[3] << ","
//	//	<< BodyBounds[4] << ","
//	//	<< BodyBounds[5] << std::endl;
//
//	BodyClipBound->SetBounds(BodyBounds);
//	vtkClipPolyData* BodyBoundsClip = vtkClipPolyData::New();
//	BodyBoundsClip->SetInputData(norma->GetOutput());
//	BodyBoundsClip->SetClipFunction(BodyClipBound);
//	BodyBoundsClip->InsideOutOn();
//	BodyBoundsClip->Update();
//	//CheckModel(BodyBoundsClip->GetOutput());
//	///
//
//
//	/////* Step 3. */
//	//vtkPolyDataNormals* NormalFilterClip = vtkPolyDataNormals::New();
//	//NormalFilterClip->SetInputData(BodyBoundsClip->GetOutput());
//	//NormalFilterClip->Update();
//	//vtkWarpVector * warpClipPoly = vtkWarpVector::New();
//	//warpClipPoly->SetInputData(NormalFilterClip->GetOutput());
//	//warpClipPoly->SetInputArrayToProcess(0, 0, 0,
//	//	vtkDataObject::FIELD_ASSOCIATION_POINTS,
//	//	vtkDataSetAttributes::NORMALS);
//	//warpClipPoly->SetScaleFactor(0.1);
//	//warpClipPoly->Update();
//	//NormalFilterClip->SetInputData(warpClipPoly->GetOutput());
//	//NormalFilterClip->Update();
//	////CheckModel(warpClipPoly->GetPolyDataOutput());
//	/////
//
//	///* Step 4. */
//	vtkBooleanOperationPolyDataFilter* booleanStep1 =
//		vtkBooleanOperationPolyDataFilter::New();
//	booleanStep1->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);
//	booleanStep1->SetInputData(0, NormalFilter->GetOutput());
//	booleanStep1->SetInputData(1, BodyBoundsClip->GetOutput());
//	booleanStep1->Update();
//
//	//TmpPolydata->DeepCopy(booleanStep1->GetOutput());
//	///
//
//	///////////////////////////////////////****************************************** HOLE Clipping
//	///* Step 5. */
//	double HoleBounds[6];
//	mpHolePolydata->GetBounds(HoleBounds);
//
//	for (int i = 0; i < 6; i++)
//	{
//		if (HoleBounds[i])
//			HoleBounds[i] *= 1.5;
//		else
//			HoleBounds[i] *= -1.5;
//	}
//
//	vtkCubeSource* cubeSource2 = vtkCubeSource::New();
//	cubeSource2->SetBounds(HoleBounds);
//	cubeSource2->SetCenter(mpHolePolydata->GetCenter());
//	cubeSource2->Update();
//
//	vtkBox* HoleClipBound = vtkBox::New();
//	HoleClipBound->SetBounds(cubeSource2->GetOutput()->GetBounds());
//
//	vtkClipPolyData* HoleBoundsClip = vtkClipPolyData::New();
//	HoleBoundsClip->SetInputData(booleanStep1->GetOutput());
//	HoleBoundsClip->SetClipFunction(HoleClipBound);
//	HoleBoundsClip->InsideOutOn();
//	HoleBoundsClip->Update();
//	HoleBoundsClip->GetOutput()->GetBounds(HoleBounds);
//
//	//vtkSphereSource* HoleSource = vtkSphereSource::New();
//	//HoleSource->SetCenter(HoleBoundsClip->GetOutput()->GetCenter());
//	//HoleSource->SetRadius((mpHolePolydata->GetBounds()[1] - mpHolePolydata->GetBounds()[0]) / 1.7);
//	//HoleSource->SetPhiResolution(150);
//	//HoleSource->SetThetaResolution(150);
//	//HoleSource->Update();
//	//CheckModel(HoleBoundsClip->GetOutput());
//	///
//
//	///* Step 6. */
//	vtkBooleanOperationPolyDataFilter* booleanStep7 = vtkBooleanOperationPolyDataFilter::New();
//	booleanStep7->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);//Intersection
//	booleanStep7->SetInputData(0, HoleBoundsClip->GetOutput());
//	booleanStep7->SetInputData(1, mpHolePolydata);
//	booleanStep7->Update();
//	///
//
//	///* Step 7. */
//	HoleBoundsClip->InsideOutOff();
//	HoleBoundsClip->Update();
//
//	vtkAppendPolyData* AppendData = vtkAppendPolyData::New();
//	AppendData->AddInputData(HoleBoundsClip->GetOutput());
//	AppendData->AddInputData(booleanStep7->GetOutput());
//	AppendData->Update();
//	//CheckModel(AppendData->GetOutput());
//	//TmpPolydata->DeepCopy(AppendData->GetOutput());
//	///
//
//	///////////////////////////////////////****************************************** SLEEVE Clipping
//	///* Step 8. */
//	vtkBox* SleeveClipBound = vtkBox::New();
//	double SleeveBounds[6];
//	mpSleevePolydata->GetBounds(SleeveBounds);
//	SleeveClipBound->SetBounds(SleeveBounds);
//	vtkClipPolyData* SleeveBoundsClip = vtkClipPolyData::New();
//	SleeveBoundsClip->SetInputData(AppendData->GetOutput());
//	SleeveBoundsClip->SetClipFunction(SleeveClipBound);
//	SleeveBoundsClip->InsideOutOn();
//	SleeveBoundsClip->Update();
//	//CheckModel(SleeveBoundsClip->GetOutput());
//
//	//std::cerr << "Step =8 Done"  << std::endl;
//
//	///
//
//	///*Step 9. 변경*/
//
//	//	//vtkNew(vtkWarpVector, warpClipPoly);
//	//vtkWarpVector* warpPolyData = vtkWarpVector::New();
//	//warpPolyData->SetInputData(TmpPolydata);
//	//warpPolyData->SetInputArrayToProcess(0, 0, 0,
//	//		vtkDataObject::FIELD_ASSOCIATION_POINTS,
//	//		vtkDataSetAttributes::NORMALS);
//	//warpPolyData->SetScaleFactor(0.01);
//	//warpPolyData->Update();
//	//	std::cerr << "Step =9-1 Done" << std::endl;
//	//CheckModel(warpPolyData->GetPolyDataOutput());
//	vtkTransform* Trans1 = vtkTransform::New();
//	Trans1->Translate(0., .0, -.1);
//
//	/*Max Tansform 적용*/
//	vtkNew(vtkTransformPolyDataFilter, xfmPd);
//	xfmPd->SetInputData(confilter->GetOutput());
//	xfmPd->SetTransform(Trans1);
//	xfmPd->Update();
//
//	vtkBooleanOperationPolyDataFilter* booleanStep4 =
//		vtkBooleanOperationPolyDataFilter::New();
//	booleanStep4->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE);//Intersection
//																				  //booleanStep4->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);
//	booleanStep4->SetInputData(0, mpSleevePolydata);
//	booleanStep4->SetInputData(1, xfmPd->GetOutput());
//	booleanStep4->Update();
//	//std::cerr << "Step =9-2 Done" << std::endl;
//	//CheckModel(booleanStep4->GetOutput());
//	///
//
//	///*Step 10.*/
//	vtkBooleanOperationPolyDataFilter* booleanStep5 =
//		vtkBooleanOperationPolyDataFilter::New();
//	booleanStep5->SetOperation(vtkBooleanOperationPolyDataFilter::VTK_UNION);//Intersection
//	booleanStep5->SetInputData(0, SleeveBoundsClip->GetOutput());
//	booleanStep5->SetInputData(1, booleanStep4->GetOutput());
//	booleanStep5->Update();
//	//CheckModel(booleanStep5->GetOutput());
//	//std::cerr << "Step =10 Done" << std::endl;
//	//TmpPolydata->DeepCopy(booleanStep5->GetOutput());
//	///
//
//	///* Step 11. */
//	SleeveBoundsClip->InsideOutOff();
//	SleeveBoundsClip->Update();
//
//	vtkAppendPolyData* AppendData2 = vtkAppendPolyData::New();
//	AppendData2->AddInputConnection(SleeveBoundsClip->GetOutputPort());
//	AppendData2->AddInputData(booleanStep5->GetOutput());
//	AppendData2->Update();
//	//CheckModel(AppendData2->GetOutput());
//	//std::cerr << "Step =11 Done" << std::endl;
//	TmpPolydata->DeepCopy(AppendData2->GetOutput());
//	///
//
//	/////* Step 8. */
//	////SleeveBoundsClip->InsideOutOff();
//	////SleeveBoundsClip->Update();
//	////HoleBoundsClip->SetInputData(SleeveBoundsClip->GetOutput());
//	//HoleBoundsClip->InsideOutOff();
//	//HoleBoundsClip->Update();
//	/////
//
//	/////* Step 9. */
//	//vtkAppendPolyData* AppendData = vtkAppendPolyData::New();
//	//AppendData->AddInputData(HoleBoundsClip->GetOutput());
//	//AppendData->AddInputData(booleanStep7->GetOutput());
//	////AppendData->AddInputData(booleanStep5->GetOutput());
//	//AppendData->Update();
//
//	//TmpPolydata->DeepCopy(AppendData->GetOutput());
//	///
//	//
//	//{
//	mpSG_Body = vtkSmartPointer<vtkPolyData>::New();
//	mpSG_Body->DeepCopy(TmpPolydata);
//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputData(TmpPolydata);
//	mapper->ScalarVisibilityOff();
//
//	/*vtkSmartPointer<vtkProperty> backFaces =
//	vtkSmartPointer<vtkProperty>::New();
//	backFaces->SetSpecular(0.5);
//	backFaces->SetDiffuse(.2);
//	backFaces->SetAmbient(1.0);
//	backFaces->SetAmbientColor(colorRGB(0, 0, 250));*/
//
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetColor(1., 1., 0.);
//	//actor->GetProperty()->LightingOff();
//	//actor->SetBackfaceProperty(backFaces);
//
//	mpMainRenderer->GetRenderer()->AddActor(actor);
//	//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//
//	//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//	//		vtkSmartPointer<vtkRenderWindow>::New();
//	//	renderWindow->AddRenderer(renderer);
//
//	//	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//	//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	//	interactor->SetRenderWindow(renderWindow);
//	//	renderer->AddActor(actor);
//
//	//	renderer->SetBackground(.1, .2, .4);
//
//	//	renderWindow->SetSize(500, 250);
//
//	//	renderWindow->Render();
//	//	interactor->Start();
//	//}
//	norma->Delete();
//	implicitModeller->Delete();
//	contourFilter->Delete();
//	reverseIN->Delete();
//	NormalFilter->Delete();
//	BodyClipBound->Delete();
//	BodyBoundsClip->Delete();
//	booleanStep1->Delete();
//	cubeSource2->Delete();
//	HoleClipBound->Delete();
//	HoleBoundsClip->Delete();
//	//HoleSource->Delete();
//	booleanStep7->Delete();
//	AppendData->Delete();
//	SleeveClipBound->Delete();
//	SleeveBoundsClip->Delete();
//	booleanStep4->Delete();
//	booleanStep5->Delete();
//	AppendData2->Delete();
//	//////////////////////////////////////////////
//
//
//	vtkSmartPointer<vtkSphereSource> sphereSource =
//		vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource->Update();
//
//	vtkSmartPointer<vtkPoints> selectionPoints =
//		vtkSmartPointer<vtkPoints>::New();
//
//	selectionPoints->InsertPoint(0, -0.16553, 0.135971, 0.451972);
//	selectionPoints->InsertPoint(1, -0.0880123, -0.134952, 0.4747);
//	selectionPoints->InsertPoint(2, 0.00292618, -0.134604, 0.482459);
//	selectionPoints->InsertPoint(3, 0.0641941, 0.067112, 0.490947);
//	selectionPoints->InsertPoint(4, 0.15577, 0.0734765, 0.469245);
//	selectionPoints->InsertPoint(5, 0.166667, -0.129217, 0.454622);
//	selectionPoints->InsertPoint(6, 0.241259, -0.123363, 0.420581);
//	selectionPoints->InsertPoint(7, 0.240334, 0.0727106, 0.432555);
//	selectionPoints->InsertPoint(8, 0.308529, 0.0844311, 0.384357);
//	selectionPoints->InsertPoint(9, 0.32672, -0.121674, 0.359187);
//	selectionPoints->InsertPoint(10, 0.380721, -0.117342, 0.302527);
//	selectionPoints->InsertPoint(11, 0.387804, 0.0455074, 0.312375);
//	selectionPoints->InsertPoint(12, 0.43943, -0.111673, 0.211707);
//	selectionPoints->InsertPoint(13, 0.470984, -0.0801913, 0.147919);
//	selectionPoints->InsertPoint(14, 0.436777, 0.0688872, 0.233021);
//	selectionPoints->InsertPoint(15, 0.44874, 0.188852, 0.109882);
//	selectionPoints->InsertPoint(16, 0.391352, 0.254285, 0.176943);
//	selectionPoints->InsertPoint(17, 0.373274, 0.154162, 0.294296);
//	selectionPoints->InsertPoint(18, 0.274659, 0.311654, 0.276609);
//	selectionPoints->InsertPoint(19, 0.206068, 0.31396, 0.329702);
//	selectionPoints->InsertPoint(20, 0.263789, 0.174982, 0.387308);
//	selectionPoints->InsertPoint(21, 0.213034, 0.175485, 0.417142);
//	selectionPoints->InsertPoint(22, 0.169113, 0.261974, 0.390286);
//	selectionPoints->InsertPoint(23, 0.102552, 0.25997, 0.414814);
//	selectionPoints->InsertPoint(24, 0.131512, 0.161254, 0.454705);
//	selectionPoints->InsertPoint(25, 0.000192443, 0.156264, 0.475307);
//	selectionPoints->InsertPoint(26, -0.0392091, 0.000251724, 0.499943);
//	selectionPoints->InsertPoint(27, -0.096161, 0.159646, 0.46438);
//
//	vtkSmartPointer<vtkSelectPolyData> loop =
//		vtkSmartPointer<vtkSelectPolyData>::New();
//	loop->SetInputConnection(sphereSource->GetOutputPort());
//	loop->SetLoop(selectionPoints);
//	loop->GenerateSelectionScalarsOn();
//	loop->SetSelectionModeToSmallestRegion(); //negative scalars inside
//
//	vtkSmartPointer<vtkClipPolyData> clip = //clips out positive region
//		vtkSmartPointer<vtkClipPolyData>::New();
//	clip->SetInputConnection(loop->GetOutputPort());
//
//	vtkSmartPointer<vtkPolyDataMapper> clipMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	clipMapper->SetInputConnection(clip->GetOutputPort());
//
//	vtkSmartPointer<vtkLODActor> clipActor =
//		vtkSmartPointer<vtkLODActor>::New();
//	clipActor->SetMapper(clipMapper);
//
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	interactor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the renderer, set the background and size
//	renderer->AddActor(clipActor);
//	renderer->SetBackground(.1, .2, .4);
//
//	renderWindow->SetSize(500, 250);
//
//	renderWindow->Render();
//	interactor->Start();
//
//
//	return EXIT_SUCCESS;
//}