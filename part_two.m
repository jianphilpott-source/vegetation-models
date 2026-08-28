function [fn_change_in_biomass]=part_two()
%FUNCTION - Calculate change in biomass on basis of NEW resource values
%Called by MAIN
%Calls to fn_change_in_grass, fn_change_in_shrub
%Returns data in 'field_species' array for change in biomass per time step
%Returns growth array to account for propagules
%--------------------------------------------------------------------------
global fieldsize randp species    
% calculated variables and arrays
global field_map field_species mid_resource B_threshold top_resource growth_rate deep_resource     
% subroutines
global fn_change_in_biomass fn_change_in_grass fn_change_in_shrub
%Initialise variables specific to subroutine       
global shrub_resource grass_resource
growth_rate=zeros(fieldsize,fieldsize,species);  %Clear growth rate array
grass_resource=zeros(fieldsize,fieldsize,randp); %Clear grasses' resource array
shrub_resource=zeros(fieldsize,fieldsize,randp); %Clear shrubs' resource array
species_1=zeros(fieldsize,fieldsize,3); species_2=zeros(fieldsize,fieldsize,3);
percent_1=zeros(fieldsize,fieldsize,3); percent_2=zeros(fieldsize,fieldsize,3);
upbound=0.000001;downbound=-0.000001; %tolerance set on threshold calculations
%-------Evaportate top layer water resource--------------------------------
%excess water has been added and moved in part one - this section removes
%some proportion of water (we need a better evaporation subroutine - not
%included by order of Tony!!)
%for column=1:fieldsize
%    for row=1:fieldsize
%        top_resource(row,column,1)=top_resource(row,column,1)*0.6;   %*0.6;
%    end
%end
%------STEP ONE - Allocate resource in each layer to each plant -----------
for column=1:fieldsize
    for row=1:fieldsize
        if (field_species(row,column,1)<0)%error check - should not be zero
            field_species(row,column,1)=0; %set break point here for debug
        end
        if (field_species(row,column,2)<0)%error check - should not be zero
            field_species(row,column,2)=0; %set break point here for debug
        end
        %allocate current time step's resource to soil layer and species
        if ((field_species(row,column,1)< 1e-10)&&(field_species(row,column,2)< 1e-10))   %if cell is unvegetated...
            for i=1:randp %... resource is divided equally between grass and shrub
                grass_resource(row,column,i)=0.5*(top_resource(row,column,i)+mid_resource(row,column,i)+deep_resource(row,column,i));
                shrub_resource(row,column,i)=0.5*(top_resource(row,column,i)+mid_resource(row,column,i)+deep_resource(row,column,i));
            end
        else   %else statement where cell is vegetated (these rules are specific to black gramma and creosotebush)      
            species_1(row,column,1)=(field_species(row,column,1)*1.44)*0.133; %split below ground biomass into layers
            species_1(row,column,2)=(field_species(row,column,1)*1.44)*0.504; %Ratio of above to below ground biomass 
            species_1(row,column,3)=(field_species(row,column,1)*1.44)*0.363; %for grass is 1:1.44  (data from ecotone)     
            species_2(row,column,1)=field_species(row,column,2)*0.067;        %for shrub is 1:1  (data from ecotone)
            species_2(row,column,2)=field_species(row,column,2)*0.32;        
            species_2(row,column,3)=field_species(row,column,2)*0.613;
            if (field_species(row,column,1)>0) %Calculate percentage of roots of each species in each soil layer)
                percent_1(row,column,1)=species_1(row,column,1)/(species_1(row,column,1)+species_2(row,column,1));
                percent_1(row,column,2)=species_1(row,column,2)/(species_1(row,column,2)+species_2(row,column,2));
                percent_1(row,column,3)=species_1(row,column,3)/(species_1(row,column,3)+species_2(row,column,3));            
            end
            if (field_species(row,column,2)>0)
                percent_2(row,column,1)=species_2(row,column,1)/(species_1(row,column,1)+species_2(row,column,1));
                percent_2(row,column,2)=species_2(row,column,2)/(species_1(row,column,2)+species_2(row,column,2));
                percent_2(row,column,3)=species_2(row,column,3)/(species_1(row,column,3)+species_2(row,column,3));           
            end        
            if (field_map(row,column)==2) %Root chanelization for shrub
                if (field_species(row,column,1)<=0)
                    ABio_prop=1; %Abio is Above ground biomass proportion
                else
                    ABio_prop=field_species(row,column,2)/(field_species(row,column,1)+field_species(row,column,2));
                end
                if (field_species(row,column,2)>0) %Simulate root channelisation process for shrubs
                    for i=1:randp %if shrub dominates cell, 56% of resource is imediately channelled to mid layer (from 'Disposition of rainwater...' by Parsons)                
                        Res_move=ABio_prop*top_resource(row,column,i);
                        top_resource(row,column,i)=top_resource(row,column,i)-Res_move;
                        mid_resource(row,column,i)=mid_resource(row,column,i)+Res_move*0.32;
                        deep_resource(row,column,i)=deep_resource(row,column,i)+Res_move*0.68;
                    end
                end
            end                
            %partition resources according to root mass - water
            grass_resource(row,column,1)=percent_1(row,column,1)*top_resource(row,column,1)+percent_1(row,column,2)*mid_resource(row,column,1)+percent_1(row,column,3)*deep_resource(row,column,1);
            shrub_resource(row,column,1)=percent_2(row,column,1)*top_resource(row,column,1)+percent_2(row,column,2)*mid_resource(row,column,1)+percent_2(row,column,3)*deep_resource(row,column,1);
           %partition resources according to root mass - nitrogen
            grass_resource(row,column,2)=percent_1(row,column,1)*top_resource(row,column,2)+percent_1(row,column,2)*mid_resource(row,column,2)+percent_1(row,column,3)*deep_resource(row,column,2);
            shrub_resource(row,column,2)=percent_2(row,column,1)*top_resource(row,column,2)+percent_2(row,column,2)*mid_resource(row,column,2)+percent_2(row,column,3)*deep_resource(row,column,2);
        end %end of if else loop where cell is vegetated
        %-------ERROR CHECKING SECTION-------------------------------------
        check=sum(species_1(row,column,:))-field_species(row,column,1)*1.44; %check that total roots = above ground biomass
        if ((check>upbound)||(check<downbound))
            disp ('part two - percent_1_error')
        end
        check=sum(species_2(row,column,:))-field_species(row,column,2);%check that total roots = above ground biomass
        if ((check>upbound)||(check<downbound))
            disp ('part two - percent_2_error')
        end
        resource_chk_1=sum(grass_resource(row,column,:))+sum(shrub_resource(row,column,:));
        resource_chk_2=sum(top_resource(row,column,:))+sum(mid_resource(row,column,:))+sum(deep_resource(row,column,:));
        resource_chk_3=resource_chk_1-resource_chk_2;
        if((resource_chk_3<downbound)||(resource_chk_3>upbound))
            disp ('resource_conservation_error in resource_chk_3, part two')
        end
    end
end

%---------STEP TWO - Change in biomass of grass----------------------------
[fn_change_in_grass]=part_two_grass();
%---------CHECK AFTER GRASS UPDATE----------------------------------------
assert(all(isfinite(field_species(:))), 'field_species contains NaN/Inf after part_two_grass');
assert(all(isfinite(grass_resource(:))), 'grass_resource contains NaN/Inf after part_two_grass');
assert(all(isfinite(mid_resource(:))), 'mid_resource contains NaN/Inf after part_two_grass');
assert(all(isfinite(deep_resource(:))), 'deep_resource contains NaN/Inf after part_two_grass');

%--------STEP THREE - Change in biomass of shrub---------------------------
[fn_change_in_shrub]=part_two_shrub();
%---------CHECK AFTER SHRUB UPDATE----------------------------------------
assert(all(isfinite(field_species(:))), 'field_species contains NaN/Inf after part_two_shrub');
assert(all(isfinite(shrub_resource(:))), 'shrub_resource contains NaN/Inf after part_two_shrub');
assert(all(isfinite(mid_resource(:))), 'mid_resource contains NaN/Inf after part_two_shrub');
assert(all(isfinite(deep_resource(:))), 'deep_resource contains NaN/Inf after part_two_shrub');

%--------STEP FOUR - Reassemble mid_resource array - -----------------------
%NOTE - new growth is added to field_species in part_one to allow handling of
%propagules in the smoosh calculations
for column=1:fieldsize
    for row=1:fieldsize        
        for i=1:randp  %remaining resource (that not used in growth or maintenance) is carried down through layers           
            mid_resource(row,column,i)=grass_resource(row,column,i);
            deep_resource(row,column,i)=shrub_resource(row,column,i);            
        end        
    end
end
%---------CHECK AFTER STEP FOUR-------------------------------------------
assert(all(isfinite(mid_resource(:))), 'mid_resource contains NaN/Inf after STEP FOUR');
assert(all(isfinite(deep_resource(:))), 'deep_resource contains NaN/Inf after STEP FOUR');

%--------STEP FIVE - Resource loss from unvegetated cells -----------------
for column=1:fieldsize %Again, a fairly arbitrary section, added to prevent resource from accumulating in empty cells
    for row=1:fieldsize %figures were chosen in a meeting with John and Tony
        if ((field_species(row,column,1)<B_threshold(1))&&(field_species(row,column,2)<B_threshold(2)))
            mid_resource(row,column,1)=mid_resource(row,column,1)*0.05;
            mid_resource(row,column,2)=mid_resource(row,column,2)*0.5;
            deep_resource(row,column,1)=deep_resource(row,column,1)*0.05;
            deep_resource(row,column,2)=deep_resource(row,column,2)*0.5;
        end
    end %end of row
%---------CHECK AFTER STEP FIVE-------------------------------------------
assert(all(isfinite(field_species(:))), 'field_species contains NaN/Inf after STEP FIVE');
assert(all(isfinite(mid_resource(:))), 'mid_resource contains NaN/Inf after STEP FIVE');
assert(all(isfinite(deep_resource(:))), 'deep_resource contains NaN/Inf after STEP FIVE');

end %end of column