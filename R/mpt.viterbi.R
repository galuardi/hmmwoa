mpt.viterbi <-
function(allpost2, L, fmat, sst, D = c(100,500), D2s = .09){
#require(EBImage)

# landmask...ADD FROM SST AND MAKE SURE THE GRID IS CLIPPED!!!!
land = make.landmask(sst)$mask

 # D is a two part (from HMM) vector for the different diffusion parameters 
 # D and D2s may be taken from kftrack/ukfsst parameter estimates...
 print(sprintf('Number of days: %i\n',dim(L[[3]][3])))
  s = D*D2s;  # I think this is  the standard deviation of the D parameter
  rrow = dim(L[[3]])[1]; # sdimensions of the output array
  ccol = dim(L[[3]])[2]
  icalc = dim(L[[3]])[3]
  numnames = 1;  # This is if we have more than one likelihood source
  
  # Define transition probabilities (convolution kernels)
unc    = sqrt(2*s[1]); 
ks = ceiling(unc*10+1); 
ks = ks + mod(ks,2) + 1;
ks1  = max(15, ks); 
kern1  = gausskern(ks1, unc)

unc    = sqrt(2*s[2]);
ks = ceiling(unc*10+1); 
ks = ks + mod(ks,2) + 1;
ks2  = max(15, ks); 
kern2  = gausskern(ks2,unc);


# % Setup output struct
mpt.maplong = matrix(L$lon,ccol,rrow,byrow=T)
mpt.maplat = matrix(L$lat,ccol,rrow,byrow=F)

dlong =  diff(L$lon)[1]  # (result.maplong[1,ccol]-result.maplong[1,1])/(ccol-1);
dlat  = diff(L$lat)[1]    #(result.maplat[rrow,1]-result.maplat[1,1])/(rrow-1);

R = mapmatrix(L$lat[1],L$lon[1],dlat, dlong);  

smatrix = allpost2
smatrix[is.nan(smatrix)] = 0
# % Setup state metric
M = smatrix[,,1]

M = log(M); #% initialise

# % Find relevant position (ones with finite probability)
subject = (M!=-Inf)*1

zro = matrix(0,rrow,ccol);
theend = icalc;

Tprevx = numeric(rrow*ccol*theend);
dim(Tprevx) = c(rrow, ccol, theend)
Tprevy=Tprevx

# Tprevy = zro;
xidx = which.min((fmat[1,8]-(L$lon))^2);
yidx = which.min((fmat[1,9]-L$lat)^2);
Tprevx[xidx,yidx,1] = xidx;  # inital lon position 
Tprevy[xidx,yidx,1] = yidx;  # inital lat position 

# Total likelihood from various sources (i.e. SST, bathymetry, light etc.)
Ltotal = smatrix*0+1
for (j in 1:numnames){   
    Ltotal = Ltotal * smatrix;
}

Ltotal[is.nan(Ltotal)] = 0
Ltotal[,,1] = smatrix[,,1]

print('Starting iterations...')
print(sprintf('Day   1 -  9...')); 

  Tx = numeric(rrow*ccol*theend);
  dim(Tx) = c(rrow, ccol, theend)
  Ty=Tx


# % Viterbi algorithm
for (j in 2:(theend-1)){
    if (!mod(j,10)){
	    print(sprintf('\n done! time = %4.3f',Sys.time())); 
        print(sprintf('Day %3.0i -%3.0i...',j,j+9)); 
		time1;
		}
  Mtemp = log(zro); #% Mtemp starts with -inf
  Ttempx = -1+zro; 
  Ttempy = Ttempx;
  
  # assess how long this thing takes to run....
  time1 = Sys.time()
  # switch td.behav(j-1)  # this is where the behavior switching is taken into account
	# case 1
	if(fmat$behav[j-1]==1){
		 ks = ks1; kern = kern1;
		 }else{
		 ks = ks2; kern = kern2
		 }
	# case 2
		# ks = ks2; kern = kern2;
# end
  
  for(xx in 1:ccol){
        for(yy in 1:rrow){
            if (as.logical(subject[yy,xx])){
			 # print(subject[yy,xx])
			 # print(xx)
			 # print(yy)

                kminlat  = 1 + max(ceiling(ks/2)-yy, 0);
                kmaxlat  =     min(ks-(yy+floor(ks/2)-rrow), ks);
                kminlong = 1 + max(ceiling(ks/2)-xx, 0);
                kmaxlong =     min(ks-(xx+floor(ks/2)-ccol), ks);
                klat = kminlat:kmaxlat; 
		        klong = kminlong:kmaxlong;
                
                mminlat  = max(yy-floor(ks/2), 1);
                mmaxlat  = min(yy+floor(ks/2), rrow);
                mminlong = max(xx-floor(ks/2), 1);
                mmaxlong = min(xx+floor(ks/2), ccol);
                mlat = mminlat:mmaxlat; 
	         	mlong = mminlong:mmaxlong;

                # % Branch matrix for current position
                B = log(Ltotal[mlat,mlong,j-1] * kern[klat,klong]);  
				
                # % Total probability of the possible tracks from (x,y)
                Msub = B + M[yy,xx]; #% sum of current state and branch metric
                
                # % Get relevant area from large array
                Mupdate  = Mtemp[mlat,mlong];
                Txupdate = Ttempx[mlat,mlong];
                Tyupdate = Ttempy[mlat,mlong];
                
                # % Find states to be updated
                update = (Mupdate<Msub);
		update[is.na(update)]=F
				
                
                # % Update in small arrays
                Mupdate[update]  = Msub[update];
                Txupdate[update] = xx;
                Tyupdate[update] = yy;
                
                # % Transfer small arrays to large arrays
                Mtemp[mlat,mlong]  = Mupdate;
                Ttempx[mlat,mlong] = Txupdate;
                Ttempy[mlat,mlong] = Tyupdate;
            }
        }
	}

	#print(Sys.time()-time1)	
   # Mtemp[Mtemp==-Inf] = 0 	
    Mtemp[land==1] = -Inf;
	
    # % Swap tracks
    subject = (Mtemp!=-Inf)*1;  # WTF?
    #subject[land==1] = 0;
	#print(sum(subject))	 
	 
    for (xx in 1:ccol){
        for (yy in 1:rrow){
            if (as.logical(subject[yy,xx])){			
                Tx[yy,xx,1:(j-1)] = Tprevx[Ttempy[yy,xx],Ttempx[yy,xx],1:(j-1)];
                Ty[yy,xx,1:(j-1)] = Tprevy[Ttempy[yy,xx],Ttempx[yy,xx],1:(j-1)];
                Tx[yy,xx,j] = xx; 
		Ty[yy,xx,j] = yy;
            }
        }
    }
   
print(Sys.time()-time1)
   # % Update the state metrics
    M = Mtemp;
    # % Store current tracks
    Tprevx = Tx; 
	Tprevy = Ty;			
}

Sys.time()-time1


M[land==1] = Inf*-1

val = max(M);
ind = which.max(M)

txy = ind2sub(c(rrow, ccol),ind);
xm = txy[2]; 
ym = txy[1]

mpt.long = Tx[ym,xm,];
mpt.long_clean = mpt.long;
mpt.lat  = Ty[ym,xm,]; 
mpt.lat_clean  = mpt.lat;

mpt.lat = mpt.lat_clean = jitter(mpt.lat)
mpt.long = mpt.long_clean = jitter(mpt.long)

mpt = pixtomap(R, mpt.long_clean, mpt.lat_clean)
mpt[,1] = mpt[,1]
mpt[1,] = as.numeric(fmat[1,8:9])
mpt[nrow(mpt),] = as.numeric(fmat[nrow(fmat),8:9])

mpt
}
