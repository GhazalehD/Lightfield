function s = Point_3Dlocation_estimate( rays )

    K = size( rays,2 );
    
    A = zeros( 2*K, 3 );
    b = zeros( 2*K, 1 );
    for k=1:K
        
        A( 2*k-1, 1 ) = rays( 1,k );
        A( 2*k-1, 2 ) = 1.0;
        b( 2*k-1 )    = rays( 3,k );
        
        A( 2*k, 1 ) = rays( 2,k );
        A( 2*k, 3 ) = 1.0;
        b( 2*k )    = rays( 4,k );

    end
    s = A \ b;
   
    
end
