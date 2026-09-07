function [fn_r_and_p]=part_one()
%FUNCTION - PART ONE - CALC CHANGE IN RESOURCE DUE TO ACTION OF VECTORS
%Called by MAIN
%Calls to fn_smoosh_water, fn_smoosh_wind, fn_smoosh_cow
%Returns data in 'Resource' array for redistributed resource per time step
%Modifies field_species array to account for propagules, represented as new growth
%--------------------------------------------------------------------------
% variables from input file (global values that are not required by the function need not be defined)
global fieldsize randp species QV 
% calculated variables and arrays
global field_species B_threshold top_resource 
global field_map growth_rate 
% functions  
global fn_r_and_p fn_smoosh_water fn_smoosh_wind fn_smoosh_cow 
%------Define variables specific to subroutine-----------------------------
upbound=0.00001;downbound=-0.00001; %upper and lower tolerance on preservation of mass check
top_resource=zeros(fieldsize,fieldsize,randp);  %Reset resource layer - resources in top 10cms are not carried over - reset here so data can be plotted in MAIN
field_map=zeros(fieldsize,fieldsize);   %clear fieldmap - map is used to identify predominant vegetation type
growth_check_1=0;growth_check_2=0;growth_check_3=0; %variables used to capture any errors
%------STEP ONE -  Map field ----------------------------------------------
 for col=1:fieldsize    %map field loop 
    for row=1:fieldsize %Map species in field - used in subsequent calcs
        [C veg]=max(field_species(row,col,:));  %find most abundant species (default is 1)
        C=any(field_species(row,col,:));   %find any zero values (C=logical=1 else 0=empty)   
        field_map(row,col)=C*veg;   %make field map with 0=empty, 1=species 1, 2=species 2
        %NOTE: new growth is added to field species HERE and not in part
        %two - long story, but in a nutshell, having it here makes it easier to capture any errors.
        for i=1:species% Growth rate handling (NB growth rate is our effective propagule term)            
            if (field_species(row,col,i)<B_threshold(i)) %if biomass is below threshold, no propagules are moved
                field_species(row,col,i)=field_species(row,col,i)+growth_rate(row,col,i);
                growth_rate(row,col,i)=0;      % set growth rate array to zero - no propagules     
            elseif (growth_rate(row,col,i)<0) %captures a negative growth rate - no propagules are moved
                field_species(row,col,i)=field_species(row,col,i)+growth_rate(row,col,i);
                growth_rate(row,col,i)=0;      % set growth rate array to zero - no propagules          
            end %end of if test for growth rate            
            if (field_species(row,col,veg)<0) %error check to capture negative biomass levels
                field_species(row,col,veg)=0; %enable breakpoint here for debugging (shouldn't ever happen)
            end %end of if test for field species             
        end %end of for loop
    end %end of row loop (map field)                                       
 end %end of col loop (map field)
 %check section to ensure conservation of biomass within calculation - find total biomass
species_check_1=sum(sum(field_species(:,:,1)))+sum(sum(growth_rate(:,:,1)))+sum(sum(field_species(:,:,2)))+sum(sum(growth_rate(:,:,2)));
%--------grazing effects--------------------------------------------------- 
% Example code for removal of grass - can be a constant value (as below) or remove in a pattern.  In either case, do it here!  
%   remove_grass=0.141;    %set value for percentage of grass to be removed
%    for col=1:fieldsize
%       for row=1:fieldsize
%           field_species(row,col,1)=field_species(row,col,1)*(1-remove_grass);
%       end
%   end
%end
%-------STEP TWO RandP Smoosh by WATER------------------------------------
[fn_smoosh_water]=part_one_water();     %Redistribution of resource and propagules by water
 %-------STEP THREE RandP Smoosh by WIND-----------------------------------
[fn_smoosh_wind]=part_one_wind();       %Redistribution of resource and propagules by wind
 %-------STEP FOUR RandP Smoosh by COWS------------------------------------
[fn_smoosh_cow]=part_one_cow();         %Redistribution of resource and propagules by cows
%-----------------STEP FIVE -  Map out smooshed resource-------------------
for i=1:randp  %combine resource parts of smoosh functions
    top_resource(:,:,i)=fn_smoosh_water(:,:,i)+fn_smoosh_wind(:,:,i)+fn_smoosh_cow(:,:,i); %this line is used when all 3 vectors operate
   %More example code lines - this is where you should control the number and type of vectors that are operating
    %top_resource(:,:,i)=fn_smoosh_water(:,:,i)+fn_smoosh_wind(:,:,i); %this line (and its variations) are used when 2 vectors operates
    %top_resource(:,:,i)=fn_smoosh_wind(:,:,i); %this line (and its variations is used when one vector operates
end
%for col=1:fieldsize %use this loop when one or more vector operates to add the rest of the resource to the top_resource array
%   for row=1:fieldsize
%        top_resource(row,col,1)=top_resource(row,col,1)+QV(1); %Assuming only wind operated, so we're adding the water resource
%        top_resource(row,col,2)=top_resource(row,col,2)+QV(2)*0.1+QV(2)*0.45;  %again, assuming wind alone, we're adding the portions that should have been moved by water and cows
%    end
%end
%-----------------STEP SIX -  Handle propagules ---------------------------
for veg=1:species %combine biomass parts of smoosh functions
    i=randp+veg;
    field_species(:,:,veg)=field_species(:,:,veg)+(fn_smoosh_water(:,:,i)+fn_smoosh_wind(:,:,i)+fn_smoosh_cow(:,:,i)); %use this line for multiple vectors
    %field_species(:,:,veg)=field_species(:,:,veg)+fn_smoosh_wind(:,:,i); %use this line if one vector operates
end %error check - conserved increase in biomass
%-----------------ERROR CHECKING (check for mass conservation)-------------
%NB - this section will return errors if the grazing section is applied,
%but shouldn't return an error if only 1 or 2 vectors are allowed to operate
tot1=zeros(fieldsize,fieldsize); tot2=zeros(fieldsize,fieldsize);
for j=1:species     %check for error in growth rate
    tot1(:,:)=(fn_smoosh_water(:,:,3)+fn_smoosh_wind(:,:,3)+fn_smoosh_cow(:,:,3));
    tot2(:,:)=(fn_smoosh_water(:,:,4)+fn_smoosh_wind(:,:,4)+fn_smoosh_cow(:,:,4));
    diff1=sum(sum(growth_rate(:,:,1)))-sum(sum(tot1(:,:)));
    diff2=sum(sum(growth_rate(:,:,2)))-sum(sum(tot2(:,:)));
    if ((diff1>upbound)||(diff1<downbound))     %I'm using bounded values as you will get numerical smearing in Matlab, albeit v small
        error_in_propagules=diff1 %enable breakpoint here for debugging
    end
    if ((diff2>upbound)||(diff2<downbound))
        error_in_propagules=diff2 %enable breakpoint here for debugging
    end
end %end of for loop  
tot3=zeros(randp,1);
for i=1:randp       %Check for mass conservation
    tot3(i)=QV(i)-(((sum(sum(top_resource(:,:,i))))/(fieldsize*fieldsize)));
    if ((tot3(i)>upbound)||(tot3(i)<downbound))
        disp ('error_in_part_one=mass_not_conserved') %enable breakpoint here for debugging
    end
end %end of for loop
species_check_2=sum(sum(field_species(:,:,1)))+sum(sum(field_species(:,:,2)));
species_check_3=species_check_1-species_check_2; %test that increase in biomass was conserved
if ((species_check_3>upbound)||(species_check_3<downbound))
    error_in_part_one_growth_rate_conservation=species_check_3
end