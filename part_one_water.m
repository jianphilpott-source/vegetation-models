function [fn_smoosh_water]=part_one_water()
%FUNCTION - smoosh_water- CALC CHANGE IN RESOURCE DUE TO ACTION OF WATER
%Called by fn_r_and_p
%Calls to fn_smoosh_resource
%Returns data in fn_smoosh_water
%--------------------------------------------------------------------------
% variables from input file
global fieldsize randp species Bmax QV Reprod 
% calculated variables and arrays
global field_species Swater B_threshold field_map growth_rate
% functions  
global fn_smoosh_water fn_smoosh_resource  
% Initialise arrays specific to subroutine
global nextrow lastrow nextcol lastcol flo_row flo_col
global rsrc  down_grad
res_and_prop=randp+species; %integer number controlling size of array (specific to this function)
rsrc=zeros(fieldsize,fieldsize,res_and_prop); %clear array - contains randp specific to this function
flo_row=1; flo_col=1; %controls execution of loop - water flow always follows field aspect
%-------STEP ONE RandP Smoosh by WATER------------------------------------
 for col=1:fieldsize
     down_grad(1:res_and_prop)=0; %randp values passed down gradient
     for row=1:fieldsize %Water always flows top to bottom
         %------------------resource availability calculations-------------
         veg_flag=field_map(row,col); %label cell-grass/shrub/empty  - dominant vegetation controls smoosh                       
         QVert(1)=QV(1);QVert(2)=0.45*QV(2); %proportion of vertical flux that the vector can operate on
         for i=1:randp
            rsrc(row,col,i)=rsrc(row,col,i)+QVert(i); %add vertical flux to cell
         end
         if (veg_flag==0||(field_species(row,col,veg_flag)<B_threshold(veg_flag))) %if cell is unvegetated
             A(1)=QVert(1); %all of the vertical flux resource is Available (A) to the vector
             A(2)=QVert(2);  
             Smoosh=zeros(3,3); %different smoosh array when cell is empty
             Smoosh(3,2)=1;  %all resource flows downslope
         else
             A(1)=(QVert(1)-(QVert(1)/Bmax(veg_flag))*field_species(row,col,veg_flag));% Availability is a function of biomass in cell
             A(2)=(QVert(2)-(QVert(2)/Bmax(veg_flag))*field_species(row,col,veg_flag));% A(i) varies inversely with Biomass                    
             Smoosh=Swater(:,:,veg_flag); 
             if (A(1)<0), A(1)=0; end %capture error when availability goes negative
             if (A(2)<0), A(2)=0; end %this error will happen when Biomass in Cell=Bmax - use for debugging
         end 
         %------------------propagule calculations-------------------------
         for veg=1:species %propagule movement is directly proportional to biomass in cell
             i=randp+veg;  %increased biomass allows an increased proportion of new growth to move to neighbour cells
             rsrc(row,col,i)=rsrc(row,col,i)+growth_rate(row,col,veg)*Reprod(veg,1);%No distinction between method e.g. diffusion (cloning, tillerage) or propagules
             A(i)=(growth_rate(row,col,veg)*Reprod(veg,1))/Bmax(veg)*field_species(row,col,veg); %vector instead accounts for a proportion of plant growth movement
         end
         %-------------------loop control calculations---------------------  
         nextrow=row+flo_row; lastrow=row-flo_row; %Water flows from top to bottom of field (row 1 to fieldsize)
         nextcol=col+flo_col; lastcol=col-flo_col;                    
         if (nextrow>fieldsize) %represents a periodic boundary condition                            
             nextrow=1;
         end         
         if (lastrow<1) %represents a periodic boundary condition  
             lastrow=fieldsize;
         end
         if (nextcol>fieldsize)  %columns are always wrapped as perpendicular to flux
             nextcol=1;
         elseif (nextcol<1)  
             nextcol=fieldsize;
         end
         if (lastcol>fieldsize)
             lastcol=1;
         elseif (lastcol<1)
             lastcol=fieldsize;
         end   
         %---------------------------smoosh calculations-------------------
         [fn_smoosh_resource]=smoosh_resource(row,col,A,Smoosh);   %call smoosh routine     
     end% end of row loop
 end% end of column loop          
fn_smoosh_water=rsrc; %return calculated data to calling function