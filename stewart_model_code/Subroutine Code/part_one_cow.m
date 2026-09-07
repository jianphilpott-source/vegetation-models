function [fn_smoosh_cow]=part_one_cow()
%FUNCTION - smoosh_cow- CALC CHANGE IN RESOURCE DUE TO ACTION OF COWS
%Called by fn_r_and_p
%Calls to fn_smoosh_resource
%Returns data in fn_smoosh_cow
%--------------------------------------------------------------------------
% variables from input file
global fieldsize randp species Bmax QV Reprod
% calculated variables and arrays
global field_species Scow B_threshold field_map growth_rate
% functions  
global fn_smoosh_cow fn_smoosh_resource
%initialise variables specific to subroutine
global nextrow lastrow nextcol lastcol flo_row flo_col
global rsrc  down_grad
res_and_prop=randp+species; %integer value to control size of array
rsrc=zeros(fieldsize,fieldsize,res_and_prop); %clear array - contains randp specific to subroutine
flo_row=-1; %identifies direction of movement thru field
flo_col=1;  %cows move uphill away from water sink
%-------STEP THREE RandP Smoosh by COWS------------------------------------
 for outer=1:fieldsize
     down_grad(1:res_and_prop)=0; %randp values passed down gradient
     for inner=1:fieldsize 
         row=fieldsize+1-inner; %assume water sink at bottom of field
         col=outer;             %so cows move upwards
         veg_flag=field_map(row,col); %label cell-grass/shrub/empty  - dominant vegetation controls smoosh  
         if (veg_flag==0), veg_flag=1;end
         QVert(1)=0; QVert(2)=0.1*QV(2); %proportion of vertical flux that the vector can operate on
         for i=1:randp
             rsrc(row,col,i)=rsrc(row,col,i)+QVert(i); %add vertical flux to cell
         end
         A(1)=0; %water isn't moved by cows
         A(2)=((QVert(2)/Bmax(veg_flag))*field_species(row,col,veg_flag));% Availability is a function of and proportional to biomass in cell
         if (A(2)<0), A(2)=0; end
         if ((sum(field_species(row,col,1:species))>sum(B_threshold(:)))) %if cell is vegetated above threshold             
             Smoosh=zeros(3,3);  %different smoosh array when cell is empty
             Smoosh(3,2)=1;      %all resource flows down-gradient
         else                                     
             Smoosh=Scow(:,:,veg_flag); %dominant vegetation type controls smoosh
         end
         %------------------propagule calculations-------------------------
         for veg=1:species
             i=randp+veg;
             rsrc(row,col,i)=rsrc(row,col,i)+growth_rate(row,col,veg)*Reprod(veg,3); %No distinction between method e.g. diffusion (cloning, tillerage) or propagules
             A(i)=(growth_rate(row,col,veg)*Reprod(veg,3))/Bmax(veg)*field_species(row,col,veg); %vector instead accounts for a proportion of plant growth movement
         end
         %----------------loop control calculations---------------------
         nextrow=row+flo_row; lastrow=row-flo_row; % Cows graze from bottom to top of field
         nextcol=col+flo_col; lastcol=col-flo_col;                            
         if (nextrow<1)  
             nextrow=fieldsize;
             %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(3,:)); %boundary condition at lower (south) end of field
             %Smoosh(3,:)=0; nextrow=fieldsize;
         end
         if (lastrow>fieldsize)
             lastrow=1;
            %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(1,:)); %boundary condition at upper (north) end of field
            %Smoosh(1,:)=0; lastrow=1;
         end           
         if (nextcol>fieldsize)   %columns are always wrapped as perpendicular to flux
             nextcol=1;
         end
         if (lastcol<1)
             lastcol=fieldsize;
         end    
         %---------------------------smoosh calculations-------------------   
         [fn_smoosh_resource]=smoosh_resource(row,col,A,Smoosh);                     
     end%end of row loop
 end%end of column loop
fn_smoosh_cow=rsrc; %return calculated variables