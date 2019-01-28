function [erpct,RMSE,avgY,avgerrorY,Bfin,oth_ind,testdatay2,Ypred2]=MC_crossval(Xfdata,yfdata,MCS,No_of_calib_sampl,pre,npcf,plotpc,varargin)
% For multiple grouped data. this code can generate subset of random sample 
% with given number of samples from each group.it will plot estimated plot
% and actual values vs its sample number in data set. Also this code will
% allow to choose any specific ranges of variables from all the variables given. 
%Xfdata , yfdata are x and  matix required for PLS
%MCS is number of montecarlo simulations to be generated as per requirement,
%it'll take scalar value. eg. 100
%min= minimum, minimum size of subset used for generating PLS model selected from whole data.
%max=maximum size of subset used for generating PLS model selected from whole data. 
%pre=preprocessing method it can handel 'mean' for mean centering and
%'auto' autoscaling.
%npcf=number of principal component to be calculated.
%plotpc plot of y estimated and actual vs sample number for that particular pc e.g. 
% if 'plotpc' =2 it will give plot of estimated y obtained by using first two principal components.
%fig_pl==1 will give plots for any other value it will not generate figures
%%%%optional inputs
%varargin{1} =datavar ='rng' if you want to select specific range of
%variables from all the variables for x data. if dont want to select range
%keep it empty ie [].varargin{2}= varrang=variable range. it will take vector that
%contains all the variables to be used, from 'x' matrix, for PLS. e.g if
%total variables are 800 and if you want to select variables from 200 to
%500 your varrng will be [200:500] or similarly [200:100 500:400] as per
%requirement.
%if data is divided in various groups based on experiments or for any other
%reasons and want to use it in this manner then varargin{3}=data and it'll take value
%data='grp' if you want to seperate data in groups samplewise. if not it'll be empty string ''
%varargin{4}=grprange. it is a vector that has range of sample number for
%different groups. eg. if you have in total 33 samples and want to divide
%them in 3 groups first grp includes samples from 1 to 19 second group
%contains samples from 20 to 24 and third group contains samples from 25 to
%30 then varargin{4}=[1 19 20 24 25 30] i.e.[start1 end1 start2 end2 start3 end3].it can have n number of groups.
%varargin{5}=grpsize. it is a vector of size n-1, if groups are n, each corrosponding to
%each group. it indicates minimum number of samples we want to include from each group 
%to form final randomly generated sample of desired size. it has one to one
%corrospondence withgroup. size of last group shall be selected based on size 
%of subset i.e min and max.eg from first group we want 10 samples and from second 
%group we want 3 samples then varargin{5}=[10 3]; and size for third will
%be decided automatically.
%%%%%%Outputs
%erpct= average absolute prediction percentage error of all 100 MCs for all
%the samples.
%RMSE= Average root mean squared error for all 100 MCs for all the samples.
%avgY= Average predicted value of each sample if first "plotpc" number of
%principal components are included in model.
%avgerrorY= Average absolute prediction percentage error  for each sample
%if first "plotpc" number of principal components are included in model.

if plotpc>1
    fig_pl=1;
else
    fig_pl=0;
    plotpc=1;
end

if nargin>7&&nargin<9
    error('varrng missing provide range of variables to be selected');
elseif nargin>7
    if isempty(varargin{1})||isempty(varargin{2})
   datavar='';
   varrng=1:size(Xfdata,2);
    else
    datavar=varargin{1};
    varrng=varargin{2};
    end
end
if (nargin>7&&nargin<12)
    error('either grprange or grpsize missing, provide range of samples for each group and there sizes'); 
elseif nargin>9
    data=varargin{3};
    grprange=varargin{4};     %group ranges i.e. number of samples in eachgroup give actual number
    grpsize=varargin{5};      % grp size number of samples from each group to be included for subset building.
     
end

%clear all
%MCS=55;                 %%%number of montecarlosimulations    
%data='grp';             %%%type of data grouped  'grp'or normal '1'
if nargin<8
    datavar='';
end
if nargin<10
    data='nogrp';
    display('no grouping of data');
end
%%%check for variable range provided

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotpc=2;               %%% pc for which y needs to be plotted.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch datavar
    case 'rng'
   % varrng=127:295 ;    %%%variable range to be selected vector.
    x=Xfdata(:,varrng);
    y=yfdata;
    otherwise
        x=Xfdata;
        y=yfdata;
        display('all the variables have been selected from original data');
end
 
[mx,~]=size(x);
 ny=size(y,2);
    if strcmp(data,'grp') 
   [~,mrng]=size(grprange);    
   grpsize=[grpsize 0];
   allgrp=[(1:round(mrng/2))' grprange(1:2:mrng)' grprange(2:2:mrng)' grpsize(1:(mrng/2))' ];   %%%% all the group ranges
       
       for e=1:size(allgrp,1)
           xg{e}=x(allgrp(e,2):allgrp(e,3),:);
       end
   end
       
      its=factorial(mx)/(factorial(No_of_calib_sampl)*factorial(mx-No_of_calib_sampl));  %% max number of combination formula
        if its>=MCS || isnan(its) || its==Inf
            its=MCS;
        end
        if nargin>10
        grpsize1{size(allgrp,1)}= No_of_calib_sampl;
        for ee1=1:(size(allgrp,1)-1)
              grpsize1{ee1}=allgrp(ee1,4);
              grpsize1{size(allgrp,1)}=grpsize1{size(allgrp,1)}-grpsize1{ee1};
        end
        end
       
        ERpct = zeros(its*ny,npcf);
        ERpctmed=zeros(its*ny,npcf);
        rsme1 =zeros(its*ny,npcf);
        rsq1=zeros(its*ny,npcf);
        RPD1=zeros(its*ny,npcf);
        press1=zeros(its*ny,npcf);
      for i=1:its                                       %%its= number of iterations and it also gives 'its' number 
        rng(i^2)                                                 %%% of different combination of samples for given size. 
        if strcmp(data,'grp')
           for e1=1:size(allgrp,1)
            Subgrpsamp{e1}=randperm(size(xg{e1},1),grpsize1{e1});
            Subgrpsamp{e1}=changem1(Subgrpsamp{e1},[allgrp(e1,2):allgrp(e1,3)],[1:size(xg{e1},1)]);          %%assigns original number to selected samples
           end
            subsamp1=cell2mat(Subgrpsamp)';       %final grouping of sample
        else
            subsamp1=randperm(mx,No_of_calib_sampl)';
        end
        
        subsamp1=sort(subsamp1);
        testsamp1=delsamps([1:mx]',subsamp1);
        datax1=x(subsamp1,:);
        testdatax1=x(testsamp1,:);
        datay1=y(subsamp1,:);
        testdatay1=y(testsamp1,:);
        switch pre
            case 'no'
                [datax2]=(datax1);
            [datay2]=(datay1);
        testdatax2=(testdatax1);
            case 'mean'
        [datax2,mnsx]=mncn(datax1);
        [datay2,mnsy]=mncn(datay1);
        testdatax2=scale(testdatax1,mnsx);
            case 'auto'
        [datax2,mnsx,sigx]=zscore(datax1);
        [datay2,mnsy,sigy]=zscore(datay1);
        testdatax2=scale(testdatax1,mnsx,sigx); 
            
        end
        [~,~,~,~,~,~,Bfin]=oplspls_v2(datax2,datay2,npcf);
        
       if (size(Bfin,1)/ny)<npcf
           a1=(size(Bfin,1)/ny);
       else 
           a1=npcf;
       end
           for j1= 1:a1
                 
             ypred = testdatax2*Bfin((j1-1)*ny+1:j1*ny,:)';
            switch pre
                case 'no'
                    ypred;
                case 'mean'
             ypred = rescale(ypred,mnsy);
                case 'auto'
                    ypred=rescale(ypred,mnsy,sigy);
            end
            
           ERpct((i-1)*ny+1:i*ny,j1) = ERpct((i-1)*ny+1:i*ny,j1)+(sum(abs((ypred-testdatay1)./testdatay1*100),1)/(size(ypred,1)))';
           Ypred(:,((ny*j1)-(ny-1)):(ny*j1))=(ypred);
           eacherrorori(:,((ny*j1)-(ny-1)):(ny*j1))=(abs((ypred-testdatay1)./testdatay1*100));
           rsme1((i-1)*ny+1:i*ny,j1) = rsme1((i-1)*ny+1:i*ny,j1)+sqrt((sum((ypred-testdatay1).^2))/(size(ypred,1)))';
           SSres=sum((ypred-testdatay1).^2);
           SStot=sum((ypred-mean(testdatay1)).^2);
           rsq1((i-1)*ny+1:i*ny,j1) = rsq1((i-1)*ny+1:i*ny,j1)+(1-SSres/SStot);
           RPD1((i-1)*ny+1:i*ny,j1) = RPD1((i-1)*ny+1:i*ny,j1)+(std(testdatay1)/sqrt(SSres/size(testdatay1,1)));
           press1((i-1)*ny+1:i*ny,j1) = press1((i-1)*ny+1:i*ny,j1)+ SSres;
           
          end
        Ypred1{i}=Ypred(:,(ny*plotpc-(ny-1)):ny*plotpc);
        Ypred2{i}=Ypred;
        eacherrorori1(testsamp1(:,1),(ny*i-(ny-1)):ny*i)= eacherrorori(:,(ny*plotpc-(ny-1)):ny*plotpc);     %%this will only store prediction error after including
        testdatay2{i}=testdatay1;                                                                                                    %%%'plot pc' number of components in the model
        testsamp2{i}=testsamp1;
      end
      
      
      for j2=1:ny
        cumERpct1(j2,:) = sum(ERpct(j2:ny:its*ny,:),1)/(size(ERpct,1)/ny)  ;
        cumrmse(j2,:)=sum(rsme1(j2:ny:its*ny,:),1)/(size(rsme1,1)/ny)  ;
        cumrsq(j2,:)=sum(rsq1(j2:ny:its*ny,:),1)/(size(rsq1,1)/ny)  ;
        cumRPD(j2,:)=sum(RPD1(j2:ny:its*ny,:),1)/(size(RPD1,1)/ny)  ;
        cumpress(j2,:)=sum(press1(j2:ny:its*ny,:),1)/(size(press1,1)/ny)  ;
   end
        erpct=reshape(cumERpct1,1,ny*npcf); 
        RMSE=reshape(cumrmse,1,ny*npcf);
        Rsq=reshape(cumrsq,1,ny*npcf);
        RPD=reshape(cumRPD,1,ny*npcf);
        PRESS=reshape(cumpress,1,ny*npcf);
        oth_ind.Rsq=Rsq;
        oth_ind.RPD=RPD;
        oth_ind.PRESS=PRESS;
 %clear subsamp datax datay Ypred testsamp1 testdatay1 %eacherrorori1
 

 
 for eej=1:ny
     eacherrorY{eej}=eacherrorori1(:,eej:ny:end);
       for eei=1:its
     Listy{eej}(testsamp2{eei},eei)=Ypred1{eei}(:,eej);
   end
   for eek=1:size(Listy{eej},1)
     Listy1{eek,eej}=Listy{eej}(eek,logical(Listy{eej}(eek,:)));
     [~,avgY{eej}(eek,1),sd{eej}(eek,1)]=zscore( Listy1{eek,eej});
     eacherrorYf{eek,eej}=eacherrorY{eej}(eek,logical(eacherrorY{eej}(eek,:)));
     [~,avgerrorY{eej}(eek,1)]=zscore(eacherrorYf{eek,eej});
   end
   end
   
   
   if fig_pl==1
   %%%%figures
   for ja=1:ny
   figure

      for jc=1:its
            %plot(testsampf{j,1}{1,jc}(:,1),Ypredf{j,1}{1,jc}(:,ja),'or'); hold on;
           %plot(testsampf{j,1}{1,jc}(:,1),testdatayf{j,1}{1,jc}(:,ja),'*b');
      
      end
      hold on;
      errorbar(avgY{ja},sd{ja},'rd','Linewidth',2);
      plot(1:size(yfdata,1),yfdata,'sk','Linewidth',2);  
      str=sprintf('Y%d actual and estimated VS Sample number PCs=%d',ja,plotpc);
      title(str,'Fontsize',14.5,'FontWeight','bold');
      set(gca,'Fontsize',16);
      xlabel('Sample Number','Fontsize',15,'FontWeight','bold')
      ylabel('Actual Y (square) Estimated Y (diamond)','Fontsize',15,'FontWeight','bold')    
    end

  figure
  for k5=1:ny
             rng((k5+2)^2)
             cmap=rand(ny,3);
             plot(erpct(1,(k5:ny:ny*npcf)),'-o','Color',cmap(k5,:));
             legend('-Dynamiclegend');
             hold all;
  end
        title('PLS-Percentage Error vs PCs');
        axis([1 npcf 0 40])
        xlabel('Latent Variable')
        ylabel('Error percentage')
        
 
  figure
  for k5=1:ny
             rng((k5+2)^2)
             cmap=rand(ny,3);
             plot(RMSE(1,(k5:ny:ny*npcf)),'-o','Color',cmap(k5,:));
             legend('-Dynamiclegend');
             hold all;
  end
        title('PLS-RMSE vs PCs');
        axis([1 npcf 0 40])
        xlabel('Latent Variable')
        ylabel('Error percentage')
   end
 