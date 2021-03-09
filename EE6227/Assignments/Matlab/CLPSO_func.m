% function [gbest_position, gbest_fitness, fitcount] = CLPSO_func(fhd, Dimension, Particle_Number, Max_Gen, Max_FES, VRmin, VRmax, Inertia_Weight, C_set, varargin)
% 	% Predefine some variables in PSO
% 	me  = Max_Gen;
% 	ps  = Particle_Number;
% 	D   = Dimension;
% 	iwt = Inertia_Weight;
% 	cc  = C_set;
	
% 	% Predefine some variables in CLPSO
% 	% Probability of Pc. Pc belongs to [0.05, 0.5]
% 	for i = 1:ps
%     	Pc(i) = 0.05 + 0.45 * (exp(10*(i-1)/(ps-1))-1)/(exp(10)-1);
% 	end
% 	% Refreshing gap
% 	m = 7;
% 	mm = 0.*ones(ps,1);

% 	fitcount = ps;
% %%-------------- initialization -----------------------

% 	mv    = 0.2 * (VRmax - VRmin);
% 	VRmin = repmat(VRmin, ps, D);
% 	VRmax = repmat(VRmax, ps, D);
% 	Vmin  = repmat(-mv, ps, D);
% 	Vmax  = -Vmin;

% 	% initialize the position of the particles
% 	pos_set = VRmin + (VRmax - VRmin) .* rand(ps, D);

% 	% initialize the velocity of the particles
% 	vel_set = Vmin  + (Vmax  - Vmin)  .* rand(ps, D);

% 	% initialize the pbest and the pbest's fitness value
% 	pbest_position = pos_set;
% 	e = feval(fhd, pos_set', varargin{:});
% 	% for i = 1:ps
% 	% 	e(i, 1) = feval(fhd, pos_set(i, :), varargin{:});
% 	% end
% 	pbest_fitness = e;

% 	% initialize the pbest, the pbest's fitness value and the gbest_position_rep
% 	[gbest_fitness, gbest_index] = min(pbest_fitness);
% 	gbest_position = pbest_position(gbest_index, :);
% 	gbest_position_rep = repmat(gbest_position, ps, 1);

% 	flag = zeros(ps, 1);
% 	% During the iteration process, the number of times that the particle doesn't change
% 	stay_num = zeros(ps, 1); 

% 	ai = zeros(ps, D);
% 	f_pbest = 1:ps;
% 	f_pbest = repmat(f_pbest', 1, D);
% 	for i = 1:ps
% 		ar = randperm(D);
% 		ai(i, ar(1:mm(i))) = 1;
% 		fi1 = ceil(ps * rand(1, D));
% 		fi2 = ceil(ps * rand(1, D));

% 		fi  = (pbest_fitness(fi1) < pbest_fitness(fi2)) .* fi1 + ... 
% 			  (pbest_fitness(fi1) >= pbest_fitness(fi2)) .* fi2;
% 		bi = ceil(rand(1, D) - 1 + Pc(i));
% 		if bi == zeros(1, D);
% 			rc = randperm(D);
% 			bi(rc(1)) = 1;
% 		end
% 		f_pbest(i, :) = bi .* fi + (1 - bi) .* f_pbest(i, :);

% 	end
% 	%上面这一段怎么感觉是 > Pc则选择fi的Pbest某一维。

% %%-------------- End of initialization -----------------------

% %iteration counter
% iter = 1;

% while iter <= me & fitcount <= Max_FES
% 	iter = iter + 1;

% 	for i = 1:ps

% 		% if mod(flag(i), m) == 0;
% 		if stay_num(i) >= m;
% 			stay_num(i) = 0;
% 			ai(i,:)=zeros(1,D);
% 			f_pbest(i, :) = i .* ones(1, D)

% 			ar = randperm(D);
% 			ai(i, ar(1:mm(i))) = 1;
% 			fi1 = ceil(ps * rand(1, D));
% 			fi2 = ceil(ps * rand(1, D));
% 			fi  = (pbest_fitness(fi1) < pbest_fitness(fi2))' .* fi1 + ... 
% 				  (pbest_fitness(fi1) >= pbest_fitness(fi2))' .* fi2;
% 			bi = ceil(rand(1, D) - 1 + Pc(i));
% 			if bi == zeros(1, D);
% 				rc = randperm(D);
% 				bi(rc(1)) = 1;
% 			end
% 			f_pbest(i, :) = bi .* fi + (1 - bi) .* f_pbest(i, :);
% 		end

% 		for dimcnt = 1:D
%             pbest_f(i, dimcnt) = pbest_position(f_pbest(i, dimcnt), dimcnt); %Update the position of the dimension Dimcnt of the particle i.
%         end

%         aa(i, :) = cc(1) .* (1 - ai(i, :)) .* rand(1, D) .* (pbest_f(i, :)) + cc(2) .* ai(i, :) .* rand(1, D) .* (gbest_position_rep(i, :) - pos_set(i, :));
        
%         vel_set(i, :) = iwt(i) .* vel_set(i, :) + aa(i, :);

%         % Check whether the velocity is within the specified range.[-mv, mv]
%         vel_set(i, :) = (vel_set(i, :) > mv) .* mv + (vel_set(i, :) <= mv) .* vel_set(i, :);
%         vel_set(i, :) = (vel_set(i, :) < (-mv)) .* (-mv) + (vel_set(i, :) >= (-mv)) .* vel_set(i, :);
% 		% Update the position of the particle i.
% 		pos_set(i, :) = pos_set(i, :) + vel_set(i, :);
% 		% pos_set(i, :)
% 		% Check if the X is within the range.
% 		% pos_set(i, :) = min(VRmax(i, :), max(pos_set(i, :), VRmin(i, :)))
% 		if (sum(pos_set(i, :) > VRmax(i, :)) + sum(pos_set(i, :) < VRmin(i, :))) == 0;
% 			e = feval(fhd, pos_set', varargin{:});
% 			fitcount = fitcount + 1;
% 			tmp = pbest_fitness(i) < e(i);
% 			if tmp == 1;
% 				stay_num(i) = stay_num(i) + 1;
% 			end
% 			temp = repmat(tmp,1,D);
% 		    pbest_position(i,:) = temp .* pbest_position(i, :) + (1-temp).*pos(i,:);
% 		    pbest_fitness(i) = tmp .* pbest_fitness(i) + (1 - tmp) .* e(i); %update the pbest
% 		    if pbest_fitness(i) < gbest_fitness
% 			    gbest_position = pbest_position(i,:);
% 			    gbest_fitness = pbest_fitness(i);
% 			    gbest_position_rep = repmat(gbest_position, ps, 1);%update the gbest
%     		end
% 		end
%     end
    
%     if fitcount >= Max_FES
%     	break;
%     end

%     if (iter == me) & (fitcount < Max_FES)
%     	iter = iter - 1;
%     end
% end

function [gbest,gbestval,fitcount, gbestval_con]= CLPSO_new_func(fhd, Dimension, Particle_Number, Max_Gen, Max_FES, VRmin, VRmax, Inertia_Weight, C_set, varargin)
%[gbest,gbestval,fitcount]= CLPSO_new_func('f8',3500,200000,30,30,-5.12,5.12)
rand('state',sum(100*clock));
% me=Max_Gen;
% ps=Particle_Number;
% D=Dimension;
% cc=[1 1];   %acceleration constants
	me  = Max_Gen;
	ps  = Particle_Number;
	D   = Dimension;
	iwt = Inertia_Weight;
	cc  = C_set;

t=0:1/(ps-1):1;t=5.*t;

Pc=0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(ps))-exp(t(1)));
% Pc=0.5.*ones(1,ps);
m=0.*ones(ps,1);
% iwt=0.9-(1:me)*(0.7/me);
% iwt=0.729-(1:me)*(0.0/me);
% cc=[1.49445 1.49445]; 
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.2*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;

pos=VRmin+(VRmax-VRmin).*rand(ps,D);

e=feval(fhd,pos',varargin{:});
% for i=1:ps;
% e(i,1)=feval(fhd,pos(i,:),varargin{:});
% end

fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);
gbestval_con(1) = gbestval;

stay_num=zeros(ps,1); 

    ai=zeros(ps,D);
    f_pbest=1:ps;f_pbest=repmat(f_pbest',1,D);
    for k=1:ps       
    ar=randperm(D);
    ai(k,ar(1:m(k)))=1;
    fi1=ceil(ps*rand(1,D));
    fi2=ceil(ps*rand(1,D));
    fi=(pbestval(fi1)<pbestval(fi2)).*fi1+(pbestval(fi1)>=pbestval(fi2)).*fi2;
    bi=ceil(rand(1,D)-1+Pc(k));
    if bi==zeros(1,D),rc=randperm(D);bi(rc(1))=1;end
    f_pbest(k,:)=bi.*fi+(1-bi).*f_pbest(k,:);

    end

    stop_num=0;
    i=1;


 while i<=me&fitcount<=Max_FES
     i=i+1;
    for k=1:ps
        
    if stay_num(k)>=5
%     if round(i/10)==i/10%|stay_num(k)>=5
    stay_num(k)=0;
    ai(k,:)=zeros(1,D);
    f_pbest(k,:)=k.*ones(1,D); 
    ar=randperm(D);
    ai(k,ar(1:m(k)))=1;
    fi1=ceil(ps*rand(1,D));
    fi2=ceil(ps*rand(1,D));
    fi=(pbestval(fi1)<pbestval(fi2)).*fi1+(pbestval(fi1)>=pbestval(fi2)).*fi2;
    bi=ceil(rand(1,D)-1+Pc(k));
    if bi==zeros(1,D),rc=randperm(D);bi(rc(1))=1;end
    f_pbest(k,:)=bi.*fi+(1-bi).*f_pbest(k,:);
    end
    
    for dimcnt=1:D
        pbest_f(k,dimcnt)=pbest(f_pbest(k,dimcnt),dimcnt);
    end
    aa(k,:)=cc(1).*(1-ai(k,:)).*rand(1,D).*(pbest_f(k,:)-pos(k,:))+cc(2).*ai(k,:).*rand(1,D).*(gbestrep(k,:)-pos(k,:));%~~~~~~~~~~~~~~~~~~~~~~  
    vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); 
    vel(k,:)=(vel(k,:)>mv).*mv+(vel(k,:)<=mv).*vel(k,:); 
    vel(k,:)=(vel(k,:)<(-mv)).*(-mv)+(vel(k,:)>=(-mv)).*vel(k,:);
    pos(k,:)=pos(k,:)+vel(k,:); 
    
    if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;
    e(k)=feval(fhd,pos(k,:)',varargin{:});
    fitcount=fitcount+1;
    tmp=(pbestval(k)<=e(k));
    if tmp==1
        stay_num(k)=stay_num(k)+1;
    end
    temp=repmat(tmp,1,D);
    pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
    pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
    if pbestval(k)<gbestval
    gbest=pbest(k,:);
    gbestval=pbestval(k);
    gbestrep=repmat(gbest,ps,1);%update the gbest
    end
    gbestval_con(i) = gbestval;
    end
    
    end

% if round(i/100)==i/100
%     plot(pos(:,D-1),pos(:,D),'b*');hold on;
%     for k=1:floor(D/2)
%         plot(gbest(:,2*k-1),gbest(:,2*k),'r*');
%     end
%     hold off
%     title(['PSO: ',num2str(i),' generations, Gbestval=',num2str(gbestval)]);  
%     axis([VRmin(1,D-1),VRmax(1,D-1),VRmin(1,D),VRmax(1,D)])
%     drawnow
% end

if fitcount>=Max_FES
    break;
end
if (i==me)&(fitcount<Max_FES)
    i=i-1;
end
end
gbestval




























