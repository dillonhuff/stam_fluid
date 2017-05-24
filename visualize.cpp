#include "visualize.h"

#include "utils.h"

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageMapper.h>
#include <vtkActor2D.h>
#include <vtkImageSlice.h>

static void CreateColorImage(const int, vtkImageData*, float*);
 
// int main(int, char *[]) {
//   vtkSmartPointer<vtkImageData> colorImage = vtkSmartPointer<vtkImageData>::New();
//   CreateColorImage(colorImage);
 
//   vtkSmartPointer<vtkImageMapper> imageMapper = vtkSmartPointer<vtkImageMapper>::New();
// #if VTK_MAJOR_VERSION <= 5
//   imageMapper->SetInputConnection(colorImage->GetProducerPort());
// #else
//   imageMapper->SetInputData(colorImage);
// #endif
//   imageMapper->SetColorWindow(255);
//   imageMapper->SetColorLevel(127.5);
 
//   vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
//   imageActor->SetMapper(imageMapper);
//   imageActor->SetPosition(20, 20);
 
//   // Setup renderers
//   vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
 
//   // Setup render window
//   vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
 
//   renderWindow->AddRenderer(renderer);
 
//   // Setup render window interactor
//   vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
 
//   vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
 
//   renderWindowInteractor->SetInteractorStyle(style);
 
//   // Render and start interaction
//   renderWindowInteractor->SetRenderWindow(renderWindow);
 
//   //renderer->AddViewProp(imageActor);
//   renderer->AddActor2D(imageActor);
 
//   renderWindow->Render();
//   renderWindowInteractor->Start();
 
//   return EXIT_SUCCESS;
// }
 
void CreateColorImage(const int N, vtkImageData* image, float* densities) {
  unsigned int dim = N;
 
  image->SetDimensions(dim, dim, 1);
  image->AllocateScalars(VTK_UNSIGNED_CHAR,3);

  for(unsigned int x = 0; x < dim; x++) {
    for(unsigned int y = 0; y < dim; y++) {
	unsigned char* pixel =
	  static_cast<unsigned char*>(image->GetScalarPointer(x,y,0));
	pixel[0] = 0;
	pixel[1] = densities[IX(x, y)];
	pixel[2] = 0;
	// if(x < dim/2)
	// 	{
	// 	pixel[0] = 255;
	// 	pixel[1] = 0;
	// 	}
	// else
	// 	{
	// 	pixel[0] = 0;
	// 	pixel[1] = 255;
	// 	}
 
	// pixel[2] = 0;
      }
    }
 
  image->Modified();
}

void visualize_density(const int N, float* densities) {
  vtkSmartPointer<vtkImageData> colorImage = vtkSmartPointer<vtkImageData>::New();
  CreateColorImage(N, colorImage, densities);
 
  vtkSmartPointer<vtkImageMapper> imageMapper =
    vtkSmartPointer<vtkImageMapper>::New();
  imageMapper->SetInputData(colorImage);

  imageMapper->SetColorWindow(255);
  imageMapper->SetColorLevel(127.5);
 
  vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
  imageActor->SetMapper(imageMapper);
  imageActor->SetPosition(0, 0);
 
  // Setup renderers
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
 
  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
 
  renderWindow->AddRenderer(renderer);
 
  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
 
  vtkSmartPointer<vtkInteractorStyleImage> style =
    vtkSmartPointer<vtkInteractorStyleImage>::New();
 
  renderWindowInteractor->SetInteractorStyle(style);
 
  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  //renderer->AddViewProp(imageActor);
  renderer->AddActor2D(imageActor);

  int imageSize[3];
  colorImage->GetDimensions(imageSize);
  renderWindow->SetSize(imageSize[0], imageSize[1]);

  renderWindow->Render();
  renderWindowInteractor->Start();
 
}

