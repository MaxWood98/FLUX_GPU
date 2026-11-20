#build arguments
buildargs = -O0 -Wall -Mbounds -Mbackslash -acc -Minfo=accel -gpu=cc89

#build settings
buildsettings = -module obj

#directories
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		io_utilities_module.o\
		flux_data.o\
		flux_geometry.o\
		flux_io.o\
		flux_mesh.o\
		flux_solve.o\
		)

#object patturn rules -> for every file in "dir"/*.f90, make the file *.o in $(OBJDIR) from it
$(OBJDIR)%.o : src/%.f90
	nvfortran $(buildsettings) $(buildargs) -c $< -o $@

$(OBJDIR)%.o : io_utilities/%.f90
	nvfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: flux_link

#linking procedure
flux_link: $(OBJS) $(addprefix $(OBJDIR), flux_main.o)
	nvfortran $(buildargs) -o flux $^ $(buildsettings) -I obj 

#clean procedure 
clean: 
	rm obj/*.mod
	rm obj/*.o 