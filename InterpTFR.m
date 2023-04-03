	
TFRrep = norm(:,5000:8076) ;

[X,Y]   = meshgrid((1:size(TFRrep,2)),(1:size(TFRrep,1))) ;
	[Xq,Yq] = meshgrid((1:size(TFRrep,2)/1000:size(TFRrep,2)),(1:size(TFRrep,1))) ;
	TFRrepInterp(:,:,1) = interp2(X, Y, TFRrep, Xq, Yq, 'spline') ;
    
    mean(TFRrepInterp,3)