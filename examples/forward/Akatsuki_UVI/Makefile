CXX       = g++-14 -std=c++2a -O3

CXXFLAGS = `pkg-config --cflags eigen3`

LDFLAGS   = -fopenmp -lnetcdf_c++4

LDLIBS    = 

TARGET    = read_result
SRC       = $(TARGET).cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f $(TARGET)