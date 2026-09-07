function [fn_change_in_shrub]=part_two_shrub()
%FUNCTION - Calculate growth rate of species 2= shrub
%Called by fn_change_in_biomass
%Returns data in 'field_species' array for change in biomass per time step
%--------------------------------------------------------------------------
global fieldsize randp Growth_lim Bmax Maintenance Efficiency fail mortality drought   
% calculated variables and arrays
global field_species growth_rate     
% subroutines
global fn_change_in_shrub
global shrub_resource 
%Initialise variables specific to subroutine       
plant=2;
%--------------------------------------------------------------------------
for column=1:fieldsize  
    for row=1:fieldsize  
        neg_flag=0; %variable used to identify a negative growth rate (loss of biomass due to insufficient resource)
        bio_mass=field_species(row,column,plant); %defined to make following calculations easier to follow
        for res=1:randp
            actual_res(res)=shrub_resource(row,column,res); %makes calculations easier to follow
            maintenance_res(res)=actual_res(res)-bio_mass*Maintenance(plant,res); %parameter is used to test resource use by current biomass
            if (maintenance_res(res)<0)  %If either resource is negative (i.e. biomass requires more resource than is available to maintain itself)
                neg_flag=1; %Insufficient resource marker
            end
        end
        if (neg_flag==1) %If insufficient resource exists
            [C lim_res]= min(maintenance_res); %find which of the two resources is limiting - min function, C isn't used
            method=lim_res; %method is a variable used to control case 'switches' (see below)
        else
            for res=1:randp
                nu_growth(res)=maintenance_res(res)/Efficiency(plant,res); %Calculate possible growth rate that resource could support
            end
            [C lim_res]=min(nu_growth); %Find which resource is most limiting to new growth
            method=3; %method is a variable used to control case 'switches' (see below)
        end
        switch lower (method)   %NB Matlab switch does not break through - i.e. if case=1 and case=2, only case=1 is executed
            case 1  %water drought - insufficient water to maintain current biomass               
                field_species(row,column,plant)=bio_mass/Maintenance(plant,1)*drought(plant); %Drought modifier reduces amount of biomass that is lost during a drought 
                shrub_resource(row,column,1)=0; %All water is used, no nirogen is used
                growth_rate(row,column,plant)=0; %No new growth, no propagules
            case 2  %nitrogen drought - insufficient nitrogen to maintain current biomass
                field_species(row,column,plant)=actual_res(2)/Maintenance(plant,2); %Calculate new biomass according to the biomass that nitrogen can support
                shrub_resource(row,column,2)=0; %All nitrogen is used
                shrub_resource(row,column,1)=actual_res(1)-field_species(row,column,plant)*Maintenance(plant,1);     
                growth_rate(row,column,plant)=0;               
            case 3  %Biomass increases - either water or nitrogen limits new growth
                max_gro=bio_mass*Growth_lim(plant); %define limiting growth rate for resource use calculation (below)
                if (nu_growth(lim_res)>max_gro)
                    nu_growth(lim_res)=max_gro; %Constrain growth rate so that it doesn't exceed plants max intrinsic growth rate of plant
                end
                growth_rate(row,column,plant)=nu_growth(lim_res); %Store all new growth in growth_rate array so that it can be smooshed in part one
                for res=1:randp
                    shrub_resource(row,column,res)=actual_res(res)-bio_mass*Maintenance(plant,res)-growth_rate(row,column,plant)*Efficiency(plant,res);
                end              
            otherwise
                disp ('error in part two shrub - resource use calculation') %Can't imagine why you'd ever see this error, but better safe than sorry
        end %End of case constuction
        growth_rate(row,column,plant)=growth_rate(row,column,plant)*(1-fail(plant)); %some proportion of the propagules fail - data from ecotone
        field_species(row,column,plant)=field_species(row,column,plant)*(1-mortality(plant)); %some proportion of plants die - data from ecotone
        if (field_species(row,column,plant)>Bmax(plant)) %If the value of biomass in current cell exceeds the maximum value of biomass of that species that can exist in a 1m^2 cell (ecotone data)
            growth_rate(row,column,plant)=growth_rate(row,column,plant)+(field_species(row,column,plant)-Bmax(plant)); %Constrain biomass to its maximum value and move excess to growth_rate to be smooshed
            field_species(row,column,plant)=Bmax(plant);
        end
        for res=1:randp
            if (shrub_resource(row,column,res)<0) %error checking section - ensure grass_resource hasn't taken a negative value
                if ((shrub_resource(row,column,res)<-0.000001)||(shrub_resource(row,column,res)>0.000001))
                    disp ('correction made in part two shrub - resource less than 0')      %enable breakpoint here for debugging - this error SHOULD NOT happen
                end
                shrub_resource(row,column,res)=0; %but, there can be very small variances in shrub resource when all resource is used
            end
        end
    end    %end of row loop
end     %end of column loop