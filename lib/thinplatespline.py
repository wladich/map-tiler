# -*- coding: utf-8 -*-
import array
import math

def array_of_size(typecode, size):
    return array.array(typecode, [0] * size)

def matrixInvert( N, input_, output ):
    
#     Receives an array of dimension NxN as input.  This is passed as a one-
#     dimensional array of N-squared size.  It produces the inverse of the
#     input matrix, returned as output, also of size N-squared.  The Gauss-
#     Jordan Elimination method is used.  (Adapted from a BASIC routine in
#     "Basic Scientific Subroutines Vol. 1", courtesy of Scott Edwards.)
#
#     Array elements 0...N-1 are for the first row, N...2N-1 are for the
#     second row, etc.
#
#     We need to have a temporary array of size N x 2N.  We'll refer to the
#     "left" and "right" halves of this array.

    tempSize = 2 * N * N
    temp = array_of_size('d', tempSize)
#     First create a double-width matrix with the input array on the left
#     and the identity matrix on the right.

    for row in xrange(N):
        for col in xrange(N):
#           Our index into the temp array is X2 because it's twice as wide
#           as the input matrix.
            temp[ 2*row*N + col ] = input_[ row*N + col ] #  left = input matrix
            temp[ 2*row*N + col + N ] = 0.0               # right = 0
        temp[ 2*row*N + row + N ] = 1.0                   #   1 on the diagonal of RHS

#     Now perform row-oriented operations to convert the left hand side
#     of temp to the identity matrix.  The inverse of input will then be
#     on the right.

    for k in xrange(N):
        if (k+1 < N): # if not on the last row
            max_ = k
            for row in xrange(k +  1, N): #  find the maximum element
                if abs( temp[row*2*N + k] ) > abs( temp[max_*2*N + k] ):
                    max_ = row
            if max_ != k: #   swap all the elements in the two rows
                for col in xrange(k, 2 * N):
                    ftemp = temp[k*2*N + col]
                    temp[k*2*N + col] = temp[max_*2*N + col]
                    temp[max_*2*N + col] = ftemp

        ftemp = temp[ k*2*N + k ]
        if ( ftemp == 0.0 ): # matrix cannot be inverted
            return False
        for col in xrange(k, 2*N):
            temp[ k*2*N + col ] /= ftemp;
        for row in xrange(N):
            if row != k:
                ftemp = temp[ row*2*N + k ]
                for col in xrange(k, 2* N):
                    temp[ row*2*N + col ] -= ftemp * temp[ k*2*N + col ]

#     Retrieve inverse from the right side of temp

    for row in xrange(0, N):
        for col in xrange(N):
            output[row*N + col] = temp[row*2*N + col + N ]

    return True

class TPS(object):
    def __init__(self, points, number_of_vars):
        nof_points = self._nof_points  = len(points)
        self._nof_vars = number_of_vars
        self._dx = None
        self._dy = None
        self.type = None
        self.rhs = [array_of_size('d', nof_points+3) for _ in xrange(number_of_vars)]
        self.coef = [array_of_size('d', nof_points+3) for _ in xrange(number_of_vars)]
        self.x = array_of_size('d', nof_points)
        self.y = array_of_size('d', nof_points)
        self.u = array_of_size('d', nof_points)
        self.index = array_of_size('L', nof_points)

        for v in xrange(number_of_vars):
            for i in xrange(nof_points):
                self.rhs[v][i+3] = points[i][2][v];
                self.x[i] = points[i][0]
                self.y[i] = points[i][1]
        self._solve()

    @staticmethod    
    def base_func( x1, y1, x2, y2 ):
        if  ( x1 == x2 ) and (y1 == y2 ):
            return 0.0;
        dist  = ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 )
        return dist * math.log(dist) 

    def _solve(self):
        x = self.x
        y = self.y;
        
        # No points at all    
        if self._nof_points < 1:
            self.type = 'VIZ_GEOREF_SPLINE_ZERO_POINTS'

        # Only one point
        elif self._nof_points == 1:
            self.type = 'VIZ_GEOREF_SPLINE_ONE_POINT'
        
        # Just 2 points - it is necessarily 1D case
        elif self._nof_points == 2:
            self._dx = x[1] - x[0]
            self._dy = y[1] - y[0]
            fact = 1.0 / ( self._dx * self._dx + self._dy * self._dy )
            self._dx *= fact
            self._dy *= fact
            self.type = 'VIZ_GEOREF_SPLINE_TWO_POINTS'
        else:
            # More than 2 points - first we have to check if it is 1D or 2D case
            xmax = x[0]
            xmin = x[0]
            ymax = y[0]
            ymin = y[0]
            sumx = 0.0
            sumy= 0.0
            sumx2 = 0.0
            sumy2 = 0.0
            sumxy = 0.0
            for p in xrange(self._nof_points):
                xx = x[p]
                yy = y[p]
                xmax = max( xmax, xx )
                xmin = min( xmin, xx )
                ymax = max( ymax, yy )
                ymin = min( ymin, yy )

                sumx  += xx;
                sumx2 += xx * xx
                sumy  += yy
                sumy2 += yy * yy
                sumxy += xx * yy
            delx = xmax - xmin
            dely = ymax - ymin

            SSxx = sumx2 - sumx * sumx / self._nof_points
            SSyy = sumy2 - sumy * sumy / self._nof_points
            SSxy = sumxy - sumx * sumy / self._nof_points

            if (delx < 0.001 * dely 
                or dely < 0.001 * delx 
                or abs ( SSxy * SSxy / ( SSxx * SSyy ) ) > 0.99 ):
                self.type = 'VIZ_GEOREF_SPLINE_ONE_DIMENSIONAL'
                unused = array_of_size('L', self._nof_points)
                self._dx = self._nof_points * sumx2 - sumx * sumx
                self._dy = self._nof_points * sumy2 - sumy * sumy
                fact = 1.0 /  math.sqrt( self._dx * self._dx + self._dy * self._dy )
                self._dx *= fact
                self._dy *= fact

                for p in xrange(self._nof_points):
                    dxp = x[p] - x[0]
                    dyp = y[p] - y[0]
                    self.u[p] = self._dx * dxp + self._dy * dyp
                    unused[p] = 1

                for p in xrange(self._nof_points):
                    min_index = -1
                    min_u = 0
                    for p1 in xrange(self._nof_points):
                        if unused[p1]:
                            if min_index < 0 or self.u[p1] < min_u:
                                min_index = p1
                                min_u = self.u[p1]
                    self.index[p] = min_index
                    unused[min_index] = 0
            else:
                self.type = 'VIZ_GEOREF_SPLINE_FULL'
                # Make the necessary memory allocations

                _nof_eqs = self._nof_points + 3
                _AA = array_of_size('d', _nof_eqs * _nof_eqs )
                _Ainv = array_of_size('d', _nof_eqs * _nof_eqs )
                
                def setA(r, c, value):
                    _AA[ _nof_eqs * r + c ] = value

                def getA(r, c):
                    return _AA[ _nof_eqs * r + c ]

                def getAinv(r,c):
                    return _Ainv[ _nof_eqs * r + c ]

                # Calc the values of the matrix A

                for c in xrange(self._nof_points):
                    setA(0,c+3, 1.0)
                    setA(1,c+3, x[c])
                    setA(2,c+3, y[c])

                    setA(c+3,0, 1.0)
                    setA(c+3,1, x[c])
                    setA(c+3,2, y[c])

                for r in xrange(self._nof_points):
                    for c in xrange(r, self._nof_points):
                        setA(r+3,c+3, self.base_func( x[r], y[r], x[c], y[c] ))
                        if r != c:
                            setA(c+3,r+3 , getA(r+3,c+3))

                # Invert the matrix
                status = matrixInvert( _nof_eqs, _AA, _Ainv )

                if not status:
                    raise ValueError('TPS calculation error: there is a problem to invert the interpolation matrix')

                # calc the coefs
                for v in xrange(self._nof_vars):
                    for r in xrange(_nof_eqs):
                        self.coef[v][r] = 0.0
                        for c in xrange(_nof_eqs):
                            self.coef[v][r] += getAinv(r,c) * self.rhs[v][c]

    
    def get_point(self, src_x, src_y):
        _nof_points = self._nof_points
        leftP = 0
        rightP = 0
        x = self.x
        y = self.y
        _dx = self._dx
        _dy = self._dy
        u = self.u
        index = self.index
        rhs = self.rhs
        coef = self.coef
        _nof_vars = self._nof_vars
        vars_ = array_of_size('d', _nof_vars)
        if self.type == 'VIZ_GEOREF_SPLINE_ZERO_POINTS':
            for v in xrange(_nof_vars):
                vars_[v] = 0.0
        elif self.type == 'VIZ_GEOREF_SPLINE_ONE_POINT':
            for v in xrange(_nof_vars):
                vars_[v] = rhs[v][3]
        elif self.type == 'VIZ_GEOREF_SPLINE_TWO_POINTS':
            fact = _dx * ( src_x - x[0] ) + _dy * ( src_y - y[0] )
            for v in xrange(_nof_vars):
                vars_[v] = ( 1 - fact ) * rhs[v][3] + fact * rhs[v][4]
        elif self.type == 'VIZ_GEOREF_SPLINE_ONE_DIMENSIONAL':
            Pu = _dx * ( src_x - x[0] ) + _dy * ( src_y - y[0] )
            if Pu <= u[index[0]]:
                leftP = index[0]
                rightP = index[1]
            elif Pu >= u[index[_nof_points-1]]:
                leftP = index[_nof_points-2]
                rightP = index[_nof_points-1]
            else:
                for r in xrange(1, _nof_points):
                    leftP = index[r-1]
                    rightP = index[r]
                    if Pu >= u[leftP] and Pu <= u[rightP]:
                        break
            fact = ( Pu - u[leftP] ) / ( u[rightP] - u[leftP] )
            for v in xrange(_nof_vars):
                vars_[v] = ( 1.0 - fact ) * rhs[v][leftP+3] + fact * rhs[v][rightP+3]
        elif self.type == 'VIZ_GEOREF_SPLINE_FULL':
            for v in xrange(_nof_vars):
                vars_[v] = coef[v][0] + coef[v][1] * src_x + coef[v][2] * src_y
            for r in xrange(_nof_points):
                tmp = self.base_func( src_x, src_y, x[r], y[r] )
                for v in xrange(_nof_vars):
                    vars_[v] += coef[v][r+3] * tmp
        else:
            raise Exception("Unknown type of TPS solver: %s" % self.type)
        return vars_.tolist();
    
if __name__ == '__main__':
    def check(points, number_of_vars, point, expected):
        tps = TPS(points, number_of_vars)
        received = tps.get_point(point[0], point[1])
        print tps.type
        assert expected == received, '%s != %s' % (expected, received)

    gcps = [[1,2, [7]], [3,-4, [8]], [5, 100, [10]]]
    check(gcps, 1, [8,8], [10.745454545454544])
    
    gcps = [[0, 0, [50, 50]], [10, 10, [100, 100]]]
    check(gcps, 2, [4,5], [72.5, 72.5])
    
    gcps = [[0, 0, [50, 50]], [10, 10, [100, 100]], [0, 10, [70, 100]]]
    check(gcps, 2, [4,5], [72.0, 75.0])
    
    check([], 1, [100,200], [0]);
    
    gcps = [[1,1, [7]], [5,5, [8]], [100, 100, [10]]]
    check(gcps, 1, [200,200], [12.105263157894736])
    
    check([[3,5,[67]]], 1, [100,200], [67]);
    
