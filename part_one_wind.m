function [fn_smoosh_wind]=part_one_wind()
%FUNCTION - smoosh_wind - CALC CHANGE IN RESOURCE DUE TO ACTION OF WIND
%Called by fn_r_and_p
%Calls to fn_smoosh_resource
%Returns data in fn_smoosh_wind
%--------------------------------------------------------------------------
% variables from input file
global fieldsize randp species windir Bmax QV Reprod
% calculated variables and arrays
global field_species Swind B_threshold field_map growth_rate
% functions  
global fn_smoosh_wind fn_smoosh_resource
% Initialise variables specific to subroutine 
global nextrow lastrow nextcol lastcol flo_row flo_col
global rsrc down_grad
res_and_prop=randp+species; %integer value to control size of array
rsrc=zeros(fieldsize,fieldsize,res_and_prop); %clear array - contains randp specific to subroutine
%------- Smoosh by WIND --------------------------------------------------- 
switch windir %Establish wind direction
    case 1      %Wind blows top to bottom of field (north to south)
        startrow=0;startcol=0;                
        inc_row=1;inc_col=1;        
    case 2      %Wind blows from bottom to top of field (south to north)
        startrow=fieldsize+1;startcol=fieldsize+1;                        
        inc_row=-1;inc_col=-1;        
    case 3      %wind blows left to right (west to east)
        startrow=fieldsize+1; startcol=0;                 
        inc_row=-1;inc_col=1;        
    case 4      %wind blows right to left (east to west)
        startrow=0; startcol=fieldsize+1;                        
        inc_row=1; inc_col=-1;        
end %end of windir switch statements
for outer=1:fieldsize
    down_grad(1:res_and_prop)=0; %contains randp values passed down-gradient
        for inner=1:fieldsize             
            if (windir==1||windir==2) %calc is row then col
                row=startrow+inc_row*inner;
                col=startcol+inc_col*outer;                 
            else %(windir==3|windir==4)calc is col then row
                row=startrow+inc_row*outer;
                col=startcol+inc_col*inner;                        
            end
            %------------------resource calculations-----------------------
            veg_flag=field_map(row,col); %label cell according to dominant vegetation 
            QVert(1)=0; A(1)=0; QVert(2)=0.45*QV(2); %proportion of vertical flux that the vector can operate on
            for i=1:randp
                rsrc(row,col,i)=rsrc(row,col,i)+QVert(i); %add vertical flux to current cell
            end                               
            if (veg_flag==0||field_species(row,col,veg_flag)<B_threshold(veg_flag))    %if cell is unvegetated                    
                 A(2)=QVert(2);
                 Smoosh=zeros(3,3); %different smoosh array when cell is empty
                 if (windir==1||windir==2) %all resource flows down-gradient
                   Smoosh(3,2)=1;
                 else %(windir==3|windir==4)calc is col then row
                   Smoosh(2,3)=1;
                 end
            else 
                 Smoosh=Swind(:,:,veg_flag); %dominant vegetation controls smoosh
                 A(2)=(QVert(2)-(QVert(2)/Bmax(veg_flag))*field_species(row,col,veg_flag));
                 if (A(2)<0), A(2)=0; end                 
                 if (windir==3||windir==4) %calc is col then row
                     Smoosh=Smoosh'; %transpose smoosh array
                 end
            end 
            %------------------propagule calculations----------------------
            for veg=1:species
                i=randp+veg; 
                rsrc(row,col,i)=rsrc(row,col,i)+growth_rate(row,col,veg)*Reprod(veg,2);%No distinction between method e.g. diffusion (cloning, tillerage) or propagules
                A(i)=(growth_rate(row,col,veg)*Reprod(veg,2))/Bmax(veg)*field_species(row,col,veg);   %vector instead accounts for a proportion of plant growth movement
            end
            %----------------loop control calculations---------------------
            flo_row=inc_row; %identifies direction of movement thru field
            flo_col=inc_col; % wind can blow in any ordinal direction   
            nextrow=row+flo_row; lastrow=row-flo_row; 
            nextcol=col+flo_col; lastcol=col-flo_col; 
            if (windir==1||windir==2)  
                if (nextrow>fieldsize) 
                    nextrow=1;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(3,:)); %boundary condition at output (south in this case)
                    %Smoosh(3,:)=0; nextrow=1;
                elseif (nextrow<1)
                    nextrow=fieldsize;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(3,:)); %boundary condition at output end (north in this case)
                    %Smoosh(3,:)=0; nextrow=fieldsize;
                end
                if (lastrow>fieldsize)                    
                    lastrow=1;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(1,:)); %boundary condition at input end (south in this case)
                    %Smoosh(1,:)=0; lastrow=1;
                elseif (lastrow<1)
                    lastrow=fieldsize;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(1,:)); %boundary condition at input end (north in this case)
                    %Smoosh(1,:)=0; lastrow=fieldsize;
                end                
                if (nextcol>fieldsize)  %columns are wrapped when wind flow is perpendicular to columns
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,3));
                    %Smoosh(:,3)=0; nextcol=fieldsize;
                    nextcol=1;
                elseif (nextcol<1)  
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,3));
                    %Smoosh(:,3)=0; nextcol=1;
                    nextcol=fieldsize;
                end
                if (lastcol>fieldsize)
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,1));
                    %Smoosh(:,1)=0; lastcol=fieldsize;
                    lastcol=1;
                elseif (lastcol<1)
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,1));
                    %Smoosh(:,1)=0; lastcol=1;
                    lastcol=fieldsize;
                end    
            else %(windir==3|windir==4) rows are wrapped when wind flow is perpendicular to row
                if (nextrow>fieldsize)  
                    nextrow=1;
                elseif (nextrow<1)  
                    nextrow=fieldsize;
                end
                if (lastrow>fieldsize)
                    lastrow=1;
                elseif (lastrow<1)
                    lastrow=fieldsize;
                end
                if (nextcol>fieldsize)  %columns are always wrapped
                    nextcol=1;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,3)); %boundary condition at output end (west)
                    %Smoosh(:,3)=0; nextcol=1;
                elseif (nextcol<1)
                    nextcol=fieldsize;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,3)); %boundary condition at output end (east)
                    %Smoosh(:,3)=0; nextcol=fieldsize;
                end
                if (lastcol>fieldsize)
                    lastcol=1;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,1)); %boundary condition at input end (east)
                    %Smoosh(:,1)=0; lastcol=1;
                elseif (lastcol<1)
                    lastcol=fieldsize;
                    %Smoosh(2,2)=Smoosh(2,2)+sum(Smoosh(:,1)); %boundary condition at input end (west)
                    %Smoosh(:,1)=0; lastcol=fieldsize;
                end                   
            end            
            %------------------------smoosh calculations-------------------
            [fn_smoosh_resource]=smoosh_resource(row,col,A,Smoosh);            
        end %end of outer loop
end %end of inner loop
fn_smoosh_wind=rsrc; %return data to calling function