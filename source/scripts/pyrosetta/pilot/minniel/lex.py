
class LexicographicalIterator :
    def __init__( self, dimsizes ) :
        self.size = len(dimsizes)
        self.dimsizes = list( dimsizes )
        self.dimprods = [1] * self.size
        for i in xrange( self.size - 1, 0, -1 ) :
            self.dimprods[ i-1 ] = self.dimprods[ i ] * self.dimsizes[ i ]
        self.search_space_size = self.dimprods[ 0 ] * self.dimsizes[ 0 ]
        self.pos = [0] * self.size
        self.at_end = False

    def increment( self ) :
        i = self.size
        while i > 0 :
            i -= 1
            self.pos[ i ] += 1
            if self.pos[ i ] == self.dimsizes[ i ] :
                self.pos[ i ] = 0
            else :
                return True
        self.at_end = True
        return False

    def upper_diagonal_increment( self ) :
        for i in xrange( self.size - 1, -1, -1 ) :
            self.pos[ i ] += 1
            if ( self.pos[ i ] == self.dimsizes[ i ]  ) :
                self.pos[ i ] = 0
            else :
                beyond_end = False
                for k in xrange(i+1,self.size) :
                    self.pos[ k ] = self.pos[ i ] + k - i
                    if self.pos[k] >= self.dimsizes[k] :
                        beyond_end = True
                        break
                if beyond_end and i == 0 :
                    for k in xrange( self.size ) :
                        self.pos[ k ] = 0
                    self.at_end = True
                    return False
                elif not beyond_end :
                    return True
        self.at_end = True
        return False


    def reset( self ) :
        for i in xrange( self.size ) : self.pos[ i ] = 0
        self.at_end = False

    def upper_diagonal_reset( self ) :
        beyond_end = False
        for i in xrange( self.size ) :
            self.pos[ i ] = i
            if i >= self.dimsizes[i] :
                beyond_end = True
        if beyond_end :
            for i in xrange( self.size ) :
                self.pos[ i ] = 0
            self.at_end = True
        else :
            self.at_end = False

    def index( self ) :
        ''' return the integer index representing the state of the lex'''
        ind = 0
        for i in xrange( self.size ) :
            ind += self.pos[ i ] * self.dimprods[ i ]
        return ind
    def set_from_index( self, ind ) :
        ''' set the state of the lex given a previously computed index'''
        for i in xrange( self.size ):
            self.pos[ i ] = ind / self.dimprods[i]
            ind = ind % self.dimprods[ i ]
