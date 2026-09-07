function [fn_smoosh_resource]=smoosh_resource(row,col,A,Smoosh)
%FUNCTION - repeated section of code that handles smoosh calculation
%variables from input file
%Called by fn_smoosh_water, fn_smoosh_wind, fn_smoosh_cow
%Returns data in rsrc
%--------------------------------------------------------------------------
global  randp species windir     
% functions  
global fn_smoosh_resource 
% Initialise arrays specific to this function
global nextrow lastrow nextcol lastcol 
global rsrc down_grad
res_and_prop=randp+species;
Qmove=zeros(res_and_prop,1);
%--------------------------------------------------------------------------
for i=1:res_and_prop
    rsrc(row,col,i)=rsrc(row,col,i)-A(i)-down_grad(i); %Give the resource that will stay behind in current cell
    Qmove(i)=A(i)+down_grad(i);
    rsrc(lastrow,lastcol,i)=rsrc(lastrow,lastcol,i)+Smoosh(1,1)*Qmove(i); %resource is moved
    rsrc(lastrow,col,i)=rsrc(lastrow,col,i)+Smoosh(1,2)*Qmove(i);         %
    rsrc(lastrow,nextcol,i)=rsrc(lastrow,nextcol,i)+Smoosh(1,3)*Qmove(i); %
    rsrc(row,lastcol,i)=rsrc(row,lastcol,i)+Smoosh(2,1)*Qmove(i);         %
    rsrc(row,col,i)=rsrc(row,col,i)+Smoosh(2,2)*Qmove(i);                 %resource that stays in cell                
    rsrc(row,nextcol,i)=rsrc(row,nextcol,i)+Smoosh(2,3)*Qmove(i);
    rsrc(nextrow,lastcol,i)=rsrc(nextrow,lastcol,i)+Smoosh(3,1)*Qmove(i); %
    rsrc(nextrow,col,i)=rsrc(nextrow,col,i)+Smoosh(3,2)*Qmove(i); 
    rsrc(nextrow,nextcol,i)=rsrc(nextrow,nextcol,i)+Smoosh(3,3)*Qmove(i); %
    down_grad(i)=Smoosh(3,2)*Qmove(i);
    if (windir==3|windir==4) %part of this is available to be smooshed in next cell
        down_grad(i)=Smoosh(2,3)*Qmove(i);
    end
end %end of randp for loop in connected vegetation     