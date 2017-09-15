function lat = lattice(arg)

	if(nargin==0)
		lat.element = cell(0);
		lat = class(lat, "lattice");
	elseif(strcmp(class(arg), "lattice"))
		lat = arg;
	elseif(iscell(arg))
		lat.element = cell(length(arg));
		for i=1:length(arg)
			lat.element{i} = arg{i};
		endfor
		lat = class(lat, "lattice");
	else
		error("Can't process input.");
	endif;

endfunction
