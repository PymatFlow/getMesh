function [O1,O2,O3,O4,O5,O6,O7,O8,O9,O10]=getMesh(input)
%Function to gather mesh information from gmesh files
%---------------------------
% DEFINE SHORT-NAME OUTPUTS:
%---------------------------
% O1=msh; Mesh file from gmesh
% O2=Mdpts; Midpoint Database
% O3=cons; Various important constants
% O4=solidpts; Nodes (vertices and midpoints) on solid domain
% O5=fluidpts; Nodes (vertices and midpoints) on fluid domain
% O6=solperimeter; Nodes (vertices and midpoints) on solid/fluid bdry
% O7=outerbdrynodes; Vertices on the outer boundary
% O8=ELsolid; Global element numbers of vertices on solid domain
% O9=ELfluid; Global element numbers of vertices on fluid domain
% O10=edgetri; List of global element #'s in solid touching bdry
%--------------------------------------------------------------------------
global searchtriangles
global solidbdrynodes
disp('Gathering mesh information...')
%Choose particular mesh to use:
if input==1
    usemesh='solidfluid.msh';
elseif input==2
    
    usemesh='solidfluid2.msh';
elseif input==3
    usemesh='solidfluid3.msh';
elseif input==4
    usemesh='solidfluid4.msh';
elseif input==5
    usemesh='solidfluid5.msh';
elseif input==6
    usemesh='solidfluid6.msh';
elseif input==7
    usemesh='solidfluid7.msh';
elseif input==-1 % <----debugging with mesh by hand
    usemesh='fixedhandmesh1.msh';
elseif input==-2 % <----debugging with mesh by hand
    usemesh='fixedhandmesh2.msh';
else
    usemesh='solidfluid.msh';
end
disp([' - Using Mesh: ',usemesh])
msh = load_gmsh2(usemesh,[1 2]); % [1,2] requests variables for lines and triangles only
M=msh.nbNod; % Number of nodes in mesh
NODECO=(1/3)*msh.POS(1:M,1:2); % NODE COordinates

msh.POS(1:M,1:2)=NODECO;
L=msh.nbTriangles; % Number of elements (triangular)
ELNODE=msh.TRIANGLES(1:L,1:4); % Global node numbers for all elements; element domain
chtype=ELNODE(1,end); % get tag for domain type; assume first type listed is fluid
k=0;
j=0;
for i=1:L % get a list of elements in each domain
    if ELNODE(i,4)~=chtype
        k=k+1;
        ELsolid(k,1:3)=ELNODE(i,1:3); %coordinates of elements in solid domain
    else
        j=j+1;
        ELfluid(j,1:3)=ELNODE(i,1:3); %coordinates of elements in fluid domain
    end
end
nELsolid=k; % number of elements in the solid domain
nELfluid=j; % number of elements in the fluid domain
% The following information can be used to subdivide the boundary into segments
%global L1
%L1=msh.nbLines; % Number of element edges on the natural boundary

%global BDRYEDGES
BDRYEDGES=msh.LINES; % EDGES (node, node, ...)
[nBdryNodes,dum]=size(BDRYEDGES); % total number of nodes on a boundary
k=1;
k2=1;
chtype=BDRYEDGES(1,end);
for i=1:nBdryNodes
    if BDRYEDGES(i,end)==chtype
        outerbdrynodes(k,1)=BDRYEDGES(i,1);
        k=k+1;
    else
        solidbdrynodes(k2,1)=BDRYEDGES(i,1);
        k2=k2+1;
    end
end
%----------------------
%Determine fluid/solid nodes
%----------------------
nfluidpts=0;
nsolidpts=0;
temp=ELfluid(:,1:3);
temp2=ELsolid(:,1:3);
for k=1:M
    check=isempty(find(k==temp));
    
    if check==0
        nfluidpts=nfluidpts+1;
        fluidpts(nfluidpts)=k;
    end
    check=isempty(find(k==temp2));
    if check==0
        nsolidpts=nsolidpts+1;
        solidpts(nsolidpts)=k;
    end
end
otherDim=nfluidpts;
nsolidnodes=nsolidpts;


% If a point is in both domains, notated as 1, fluid=2, solid=0
tol=10^(-8);
Mdpts=zeros(3*M,8); %pre-allocate space since no more than triple the nodes
for mdpts
    Mdpts(:,3)=ones(3*M,1); %assume midpoint is on boundary until proven not
    nMdpts=1; % will count the number of unique midpoints in the mesh
    for i=1:L % cycle through each triangle (can we skip solid domain midpts?)
        v1=ELNODE(i,1); v2=ELNODE(i,2); v3=ELNODE(i,3);
        mid=zeros(3,4);
        mid(1,1:3)=mdptcalc(v2,v3); %organize midpoints in same way as local basis        phi's
        mid(2,1:3)=mdptcalc(v3,v1);
        mid(3,1:3)=mdptcalc(v1,v2);
        mid(:,4)=[4;5;6];
        %Check to see if any of the above are repeats
        for j=1:3
            vec=ones(nMdpts,1);
            diff_x=abs(mid(j,1)*vec-Mdpts(1:nMdpts,1));
            diff_y=abs(mid(j,2)*vec-Mdpts(1:nMdpts,2));
            check=min(diff_x+diff_y);
            if check>tol
                Mdpts(nMdpts,1:2)=mid(j,1:2);
                Mdpts(nMdpts,4)=mid(j,4);
                Mdpts(nMdpts,5)=i;
                Mdpts(nMdpts,8)=mid(j,3);
                nMdpts=nMdpts+1;
                
            else
                [dum,place]=min(diff_x+diff_y);
                Mdpts(place,6:7)=[j+3,i];
                Mdpts(place,3)=0; %mark that your midpoint isn't on bdry
            end
        end
    end
    nMdpts=nMdpts-1;
    Mdpts=Mdpts(1:nMdpts,1:8);
    %----------------------------------------
    %Determine the fluid/solid midpoints
    %----------------------------------------
    temp=Mdpts(:,8);
    count=0;
    clear solperimeter
    for k=1:nMdpts
        if temp(k)~=0
            nfluidpts=nfluidpts+1;
            fluidpts(nfluidpts)=k+M;
        end
        if temp(k)~=2
            nsolidpts=nsolidpts+1;
            solidpts(nsolidpts)=k+M;
            if temp(k)==1
                count=count+1;
                
                solperimeter(count,1)=k+M;
            end
        end
    end
    solperimeter=[solidbdrynodes; solperimeter];
    nsolperimeter=length(solperimeter);
    %Build a list of global node numbers that intersect the solid boundary,
    %edgetri
    edgetri=zeros(0,2);
    count=0;
    chtype=ELNODE(1,end); % get tag for domain type; assume first type
    listed is fluid
    for q=1:L
        if ELNODE(q,end)~=chtype % if we are looking at an element in solid
            domain
            count=count+1;
            v1=ELNODE(q,1);
            test1=isempty(find(v1==solidbdrynodes)); %see if any of the vertices
            are on solid bdry
            v2=ELNODE(q,2);
            test2=isempty(find(v2==solidbdrynodes)); %test=0 means yes on the bdry
            v3=ELNODE(q,3);
            test3=isempty(find(v3==solidbdrynodes));
            if test1+test2+test3==1 % if two are on bdry
                
                
                
                edgetri=[edgetri;[q,count]]; %record the global element number and
                solid global element number
            end
        end
    end
    nedgetri=length(edgetri); % number of solid elements on bdry
    cons=[nELsolid,nELfluid,otherDim,nsolidpts,nfluidpts,nMdpts,nsolperimeter,nedgetri];
    searchtriangles=[Mdpts(:,5),Mdpts(:,7)]; %Global element numbers for each
    midpoint
    % DEFINE SHORT-NAME OUTPUTS:
    O1=msh;
    O2=Mdpts;
    O3=cons;
    O4=solidpts;
    O5=fluidpts;
    O6=solperimeter;
    O7=outerbdrynodes;
    O8=ELsolid;
    O9=ELfluid;
    O10=edgetri;
end

%-------------------
% Midpoint Database
%-------------------
    function [out]=mdptcalc(v1,v2)
        % Gives midpoint (x,y) coordinates between global nodes v1 and v2
        % also places a flag for the domain
        % out(3)=0 means pt in interior of solid
        % out(3)=1 means pt on solid bdry
        % out(3)=2 means pt outside solid domain not on solid bdry
        v1coord=NODECO(v1,1:2);
        v2coord=NODECO(v2,1:2);
        out(1)=.5*(v1coord(1)+v2coord(1));
        out(2)=.5*(v1coord(2)+v2coord(2));
        out(3)=0;
        
        check1=isempty(find(v1==fluidpts)); %is v1 in fluid domain, 0=yes,1=no
        check2=isempty(find(v2==fluidpts)); %is v2 in fluid domain, 0=yes,1=no
        if check1+check2==0 %if both vertices are in fluid domain
            index1=find(v1==solidbdrynodes); % where is v1 on solid bdry if at             all?
            index2=find(v2==solidbdrynodes);
            check3=isempty(index1);
            check4=isempty(index2);
            if check3+check4==0 %if both vertices are on the solid bdry
                if abs(index1-index2)==1 || abs(index1-index2)==length(solidbdrynodes)-1
                    out(3)=1;
                end
            else
                out(3)=2;
            end
        end
    end
disp('Generating midpoint database...')
disp(' ')
%Mdpts= [x-coord, y-coord, bdry?, local node #1, global elem #1, local node2,
% global element #2, fluid/solid]
% loctype is from 1-6 depending on what node is across
% midpoint could be used twice, so have 2 loctype and node2