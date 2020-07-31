function y=SeparateIceMasses(lat,long)

% SeparateIceMasses separates lat and long coordinates into different ice masses
%
% 1 = Laurentide
% 2 = Greenland
% 3 = European
% 4 = West Antarctica
% 5 = East Antarctica
% 6 = Other NH (Iceland)
% 7 = Other SH (Patagonia, New Zealand, Tasmania)

 y = zeros(size(lat));

% We define Antarctica as everything below 58 S

 antarctic = (lat<-58);

% West Antarctica is from 40 W west to 173 E
% East Antarctica is the rest

% waisbndy = [  179.69  -58 ;  179.69  -67.21 ;  159.47  -79.16 ;  162.51  -82.30 ;  219.07  -86.21 ;   323.68  -80.33 ;  325.96  -78.54 ;  325.96  -58 ; 179.69 -58];

 waisbndy = [  175 -58 ;  175  -72 ;  164.5  -76.5 ;  164.5  -78 ;  219  -82 ;   302  -78 ;  305  -58 ];

 wais = inpolygon(long,lat,waisbndy(:,1),waisbndy(:,2));
 wais = antarctic .* (wais==1);
 eais = antarctic .* (wais==0);

 y(find(wais)) = 4;
 y(find(eais)) = 5;

% Northern Hemisphere Ice is north of 35 N
% b1 and b2 demarcate the boundaries of Greenland
% Europe is in longitude [0,120] and [347,360]
% North America is in longitude [185, 310]

 nhice = (lat>35);

 grnbndy = [320 55 ; 305 60 ; 300 72 ; 285 78 ; 293 81 ; 300 82.25 ; 300 90 ; 353 90 ; 353 75 ; 335 67 ; 320 55 ]; 
 grn = inpolygon(long,lat,grnbndy(:,1),grnbndy(:,2));
 y(find(grn)) = 2;

 laur = (nhice==1) .* (grn==0) .* (long>=185) .* (long<=310);
 y(find(laur)) = 1;
 
 eur = (nhice==1) .* (grn==0) .* ((long<=120) + (long>=347));
 y(find(eur)) = 3;
 y(find((y==0).*(lat>0))) = 6;
 y(find(y==0)) = 7;

end