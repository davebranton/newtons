// LibPNG example
// A.Greensted
// http://www.labbookpages.co.uk

// Version 2.0
// With some minor corrections to Mandlebrot code (thanks to Jan-Oliver)

// Version 1.0 - Initial release

#include <stdio.h>
#include <math.h>
#include <png.h>
#include <memory>
#include <vector>
#include <functional>
#include <complex>
#include <thread>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <deque>
#include <mutex>
#include <condition_variable>

#include "argstream.h"


/* Extended Array */
struct rgb
{
    unsigned char r, g, b;
};
rgb colours[] = {
    {197,197,197},
    {76,91,97},
    {130,145,145},
    {148,155,150},
    {44,66,63},
    {76,91,97},
    {130,145,145},
    {148,155,150},
};

template<class _Tnumber, class _Titerator >
bool next_combination
(
    _Titerator const& _First
    , _Titerator const& _Last
    , _Tnumber const& _Max //!< Upper bound. Not reachable
)
{
    _Titerator _Current = _First;
    if( _Current  == _Last )
    {
        return false;
    }
    *_Current += 1;
    if( *_Current < _Max )
    {
        return true;
    }
    _Titerator _Next = _Current + 1;
    if( _Next == _Last )
    {
        return false;
    }
    if( false == next_combination( _Next, _Last, _Max - 1 ) )
    {
        return false;
    }
    *_Current = *_Next + 1; 
    return *_Current < _Max;
}

int writeImage(const char* filename, int width, int height, int *buffer, const char* title);

int main(int argc, char *argv[])
{

  argstream::argstream<char> as(argc, argv);

  //constexpr int iHeight = 4096*2;
  //constexpr int iWidth = 4096*2;

  int iHeight = 1024;
  int iWidth = 1024;


  int poly_degree = 3;
  double scale = 1.0;
  std::complex<double> offset(0,0);
  double spiral_scale = 1.0 , spiral_scale_factor = 1.0;

  std::string filename;

  double real_ellipse_factor = 1.0, imag_ellipse_factor = 1.0;

  int num_cpus = 4;
  int max_iters = 20;

  double rotation = 0.0;
  double overlap = 1.0;

  as >> argstream::copyright("RootPrism");
  as >> argstream::parameter('w', "width", iWidth, "Image Width", false);
  as >> argstream::parameter('h', "height", iHeight, "Image Height", false);
  as >> argstream::parameter('o', "file", filename, "Output filename", true);
  as >> argstream::parameter('d', "degree", poly_degree, "Degree of polynomial",false);
  as >> argstream::parameter('r', "real_ellipse_factor", real_ellipse_factor, "Flatten to ellipse along real axis", false);
  as >> argstream::parameter('i', "imag_ellipse_factor", real_ellipse_factor, "Flatten to ellipse along imaginary axis", false);
  as >> argstream::parameter('s', "scale", scale, "Scale image", false);
  as >> argstream::parameter('p', "spiral_scale", spiral_scale, "Spiral scale start", false);
  as >> argstream::parameter('f', "spiral_factor", spiral_scale_factor, "Spiral scale factor", false);
  as >> argstream::parameter('n', "n_threads", num_cpus, "Number of threads to use", false);
  as >> argstream::parameter('m', "max_iters", max_iters, "Max iterations", false);
  as >> argstream::parameter('t', "rotation", rotation, "Radians to rotate roots by", false);
  as >> argstream::parameter('v', "overlap", overlap, "Root circle overlap factor", false);
  
  
  as.defaultErrorHandling();

	// Make sure that the output filename argument has been provided
	if (filename.empty()) {
		fprintf(stderr, "Please specify output file\n");
    std::cout << as.usage();
		return 1;
	}

	// Create a test image - in this case a Mandelbrot Set fractal
	// The output is a 1D array of floats, length: width * height
	printf("Creating Image\n");
    std::vector<int> buffer(iWidth * iHeight);

    // nth roots of unity
    std::vector<std::complex<double>> roots(poly_degree);
    
    for(int i = 0; i < poly_degree; i++)
    {
        roots[i] = std::exp(std::complex<double>(0, 
        rotation + (overlap*2.0*double(i)*3.141592654) / double(poly_degree)));
        roots[i].real( roots[i].real() * real_ellipse_factor );
        roots[i].imag( roots[i].imag() * imag_ellipse_factor );
        roots[i] *= spiral_scale;
        spiral_scale = spiral_scale_factor;
    }

    // from the roots, find the polynomal coeffs
    std::vector<std::complex<double>> coeffs(poly_degree+1);
    
    for(int power = 0; power < poly_degree; power++)
    {
        std::vector<int> term_indexes;
        for(int i = poly_degree - (1 + power); i >= 0; i--)
        {
            term_indexes.emplace_back(i);
        }
        do
        {
            // multiply 
            std::complex<double> coeff(1,0);
            for(auto term_index : term_indexes)
            {   
                coeff *= -roots[term_index];
            }
            coeffs[power] += coeff;
        } while (next_combination(term_indexes.begin(), term_indexes.end(), poly_degree));
    }
    coeffs[poly_degree].real(1.0);
    int ipower = 0;
    std::cout << "Polynomial coefficients" << std::endl;
    for(auto c : coeffs)
    {
        std::cout << ipower << ":" << c << std::endl;
        ipower++;
    }

    // threads
    std::vector<std::thread> threads;
    std::deque<bool> lines_done;
    std::mutex lines_lock;
    std::condition_variable cv;
    for(int div = 0; div < num_cpus; div++)
    {
        threads.emplace_back([&, div]()
        {
            std::vector<double> distances(poly_degree);
            for(int x = div*(iWidth / num_cpus); x < (div+1)*(iWidth / num_cpus) ; x++)
            {
                int last_value = 0;
                for(int y = 0; y < iHeight ; y++)
                {
                    std::complex<double> xn( double(x - iWidth/2) / double(iWidth / 2), double(y - iHeight/2) / double(iHeight / 2) );
                    xn *= scale;
                    xn += offset;
                    int i = 0;
                    do
                    {
                        // evaluate poly at xn
                        
                        std::complex<double> poly_value{};
                        for(int power = 0; power < poly_degree + 1; power++)
                        {
                            poly_value += coeffs[power] * std::pow(xn, power);
                        }

                        // derivative
                        std::complex<double> derv_value{};
                        for(int power = 1; power < poly_degree + 1; power++)
                        {
                            derv_value += coeffs[power] * std::pow(xn, power - 1) * double(power);
                        }
                        auto delta = -(poly_value / derv_value);
                        
                        //auto delta = - (std::pow(xn, poly_degree) - std::complex<double>(1,0)) / (std::pow(xn, poly_degree - 1)*(double)poly_degree);
                        if (std::norm(delta) < 0.00000001)
                        {
                          break;
                        }
                        xn += delta;
                        i++;
                        // break out toooo many
                    } while (i < max_iters);
                    // which root we closest to?
                    for(int root = 0; root < poly_degree; root++)
                    {
                        distances[root] = std::norm(xn - roots[root]);
                    }
                    auto this_value = std::min_element(distances.begin(), distances.end()) - distances.begin();
                    if (this_value != last_value)
                    {
                        buffer[y*iWidth + x] = 3;
                    }
                    else
                    {
                        buffer[y*iWidth + x] = 4;
                    }
                    last_value = this_value;
                    //std::cerr << xn << ':' << smallest << std::endl;
                }
                lines_lock.lock();
                lines_done.push_front(true);
                cv.notify_one();
                lines_lock.unlock();
            }
          lines_lock.lock();
          lines_done.push_front(false);
          cv.notify_one();
          lines_lock.unlock();
        });
    }
    threads.emplace_back([&](){
      int done = num_cpus;
      int completed = 0;
      while(done > 0)
      {
        std::unique_lock<std::mutex> ul(lines_lock);
        cv.wait(ul, [&](){
          return !lines_done.empty();
        });
        if (lines_done.back())
        {
          completed++;
          std::cout << double(int(1000.0*(double(completed) / (double(iWidth))))) / 10.0 << "%   \r";
          std::flush(std::cout);
        }
        else
        {
          done--;
        }
        lines_done.pop_back();
      };
    });
    for(auto& t : threads)
    {
        t.join();
    }
    std::cout << std::endl;

	// Save the image to a PNG file
	// The 'title' string is stored as part of the PNG file
	printf("Saving PNG\n");
	int result = writeImage(filename.c_str(), iWidth, iHeight, buffer.data(), "This is my test image");

	return result;
}

inline void setRGB(png_byte *ptr, int val)
{
    ptr[0] = colours[val].r;
    ptr[1] = colours[val].g;
    ptr[2] = colours[val].b;
}

int writeImage(const char* filename, int width, int height, int *buffer, const char* title)
{
	// Open file for writing (binary mode)
	std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(filename, "wb"), &fclose);
	if (!fp) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
        return 1;
	}

	// Initialize write structure
	std::unique_ptr<png_struct, std::function<void(png_struct*)>> png_ptr(png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL), [](png_struct*s)
    {
        png_destroy_write_struct(&s, nullptr);
    });
	if (!png_ptr) {
        return 1;
	}

	// Initialize info structure
	std::unique_ptr<png_info, std::function<void(png_info*)>> info_ptr(png_create_info_struct(png_ptr.get()), [&](png_info*p)
    {
        png_free_data(png_ptr.get(), p, PNG_FREE_ALL, -1);
    });
	if (!info_ptr) {
        return 1;
	}

	png_init_io(png_ptr.get(), fp.get());

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr.get(), info_ptr.get(), width, height,
			8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	// Set title
	if (title != NULL) {
		png_text title_text;

		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		title_text.key = (char*)"Title";
		title_text.text = (char*)title;
		png_set_text(png_ptr.get(), info_ptr.get(), &title_text, 1);
	}

	png_write_info(png_ptr.get(), info_ptr.get());

	// Allocate memory for one row (3 bytes per pixel - RGB)
	auto row = std::vector<png_byte>(3 * width);

	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			setRGB(&(row[x*3]), buffer[y*width + x]);
		}
		png_write_row(png_ptr.get(), row.data());
	}

	// End write
	png_write_end(png_ptr.get(), NULL);

	return 0;
}