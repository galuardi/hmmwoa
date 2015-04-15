mapmatrix <-
function(x, y, dx, dy){
	M = matrix(0,4,2);
	M[1,1] = y-dy;
	M[2,1] = dy;
	M[3,2] = x-dx;
	M[4,2] = dx;
	M
}
