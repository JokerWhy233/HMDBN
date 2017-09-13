function [seeddag hiddenGraph_Ps SampleDistribution] = hmdbn_structEM(TimeSeriesData, max_fan_in)

fprintf('\n');
fprintf('###########################################################################\n');
fprintf('########   Hidden Markov induced Dynamic Bayesian Network (HMDBN)  ########\n');
fprintf('########                                            by Shijia Zhu  ########\n');
fprintf('###########################################################################\n');
fprintf('\n');

    [ns ts]=size(TimeSeriesData);
    data=zeros(2*ns,(ts-1));
    data(1:ns,1:(ts-1))=(TimeSeriesData(:,1:(ts-1)));
    data((ns+1):(2*ns),1:(ts-1))=(TimeSeriesData(:,2:(ts)));
    [N ncases] = size(data);
    
    
    hiddenGraph_Ps = cell(N,1);
    SampleDistribution = cell(N,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  initial values for transition matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if nargin < 2, max_fan_in = N; end
    ns=2*ones(1,N);
    NG=2;
    init_state_distrib=ones(1,NG)/NG;
    ot=(1/ncases)/(NG-1);
    transmat=diag((1-2*ot)*ones(1,NG),0)+ot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  posterior probability for hmdbn with two nodes 
%  that is gamma in forward and backward algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('==> Calculating initial distribution Pr(qt|X,HMDBN)...\n');
    
    reversion=0;
    allgamma=cell(N,N);
    gamma=zeros(2,ncases)+0.5;
    obslik=zeros(NG,ncases);
   
    for i=1:(N/2)
       for j=(N/2+1):N
             
             tempCPDnode=cell(1,2);
             for gi=1:NG
                 if(gi==1)
                    ps=i;
                 end
                 if(gi==2)
                    ps=[];
                 end
                 fam = [ps j];
                 tempCPDnode{gi} = hmdbn_learn_params(fam,ns, data(:,:), gamma);
             end
             for gi=1:NG
                 if(gi==1)
                    ps=i;
                 end
                 if(gi==2)
                    ps=[];
                 end

                 self_ev=data(j,:);
                 pev=data(ps,:);
                 CPD= tempCPDnode{gi};
                 obslik(gi,:) = hmdbn_prob_node(CPD,self_ev, pev);
             end
             [alpha, beta, allgamma{i,j},loglik,tempxi] = hmdbn_fwdback(init_state_distrib, transmat, obslik);   
        end  
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  initial scores for hmdbn without edges   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fprintf('==> Calculating initial BWBIC score...\n');
    
    allscore=zeros(1,N); 
    for j=(1):N
        ps=[];   
        [allscore(j) temp_gamma temp_hiddenGraph_Ps]= hmdbn_hiddenGraphs_and_BWBIC_node(j, ps, ns,  data,allgamma);
    end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   greedy climing to learn hmdbn   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('==> Reconstructing the time-evolving DBN...\n');
	
    seeddag=zeros(N,N);
    outcount=0;
    while outcount < 2
        innercount = 0;

        for i=1:(N/2)
            
            fprintf('  ==> scanning Node%d...\n',i);
            
            for j=(N/2+1):N
            
               if i==j, continue;    end;
               if (i+N/2)==j, continue;    end;
               if seeddag(i,j) == 0  % No edge i-->j, then try to add it
                   tempdag = seeddag;
                   tempdag(i,j)=1;
                   if 1 % acyclic(tempdag)
                        ps=parents(tempdag,j);   
                        [temp_score  temp_gamma temp_hiddenGraph_Ps]= hmdbn_hiddenGraphs_and_BWBIC_node(j, ps, ns,  data,allgamma);

                        if temp_score > allscore(j)  &  sum( seeddag(:,j) ) < max_fan_in
                            
                            seeddag = tempdag;
                            allscore(j)= temp_score;
                            SampleDistribution{j} = temp_gamma;
                            hiddenGraph_Ps{j} = temp_hiddenGraph_Ps;
                            innercount =innercount+1;

                            fprintf('    ( + ) add edge: Node%d-->Node%d\n',i,j-N/2);
                            
                        end;
                   end
               else  % exists edge i--j, then try reverse it or remove it
                   tempdag = seeddag;
                   tempdag(i,j) = 0; 
                   tempdag(j,i) = 1;
                   if reversion % (acyclic(tempdag) & reversion) 
                        ps=parents(tempdag,i); 
                        [temp_scorei  temp_gammai temp_hiddenGraph_Psi]= hmdbn_hiddenGraphs_and_BWBIC_node(i, ps, ns,  data,allgamma);
                        ps=parents(tempdag,j);   
                        [temp_scorej  temp_gammaj temp_hiddenGraph_Psj]= hmdbn_hiddenGraphs_and_BWBIC_node(j, ps, ns,  data,allgamma);
                       if (temp_scorei+temp_scorej) > (allscore(i)+allscore(j))
                          
                           fprintf('    ( <- ) reverse edge: Node%d-->Node%d\n',i,j-N/2);
                           
                           seeddag = tempdag;
                           innercount =innercount+1;
                           allscore(i)=temp_scorei;
                           allscore(j)=temp_scorej;
                           SampleDistribution{i} = temp_gammai;
                           SampleDistribution{j} = temp_gammaj;
                           hiddenGraph_Ps{i} = temp_hiddenGraph_Psi;
                           hiddenGraph_Ps{j} = temp_hiddenGraph_Psj;
                             
                       else
                           tempdag = seeddag;
                           tempdag(i,j) = 0;
                           ps=parents(tempdag,j);   
                           [temp_score  temp_gamma temp_hiddenGraph_Ps]= hmdbn_hiddenGraphs_and_BWBIC_node(j, ps, ns,  data,allgamma);
                           if temp_score > allscore(j)
                              
                               fprintf('    ( - ) remove edge: Node%d-->Node%d\n',i,j-N/2);
                               
                               innercount =innercount+1;
                               seeddag = tempdag;
                               allscore(j)= temp_score;
                               SampleDistribution{j} = temp_gamma;
                               hiddenGraph_Ps{j} = temp_hiddenGraph_Ps;
                           
                           end;
                       end;
                   else
                           tempdag = seeddag;
                           tempdag(i,j) = 0;
                           ps=parents(tempdag,j);   
                           [temp_score  temp_gamma temp_hiddenGraph_Ps]= hmdbn_hiddenGraphs_and_BWBIC_node(j, ps, ns,  data,allgamma);
                           
                           if temp_score > allscore(j)
                              
                               fprintf('    ( - ) remove edge: Node%d-->Node%d\n',i,j-N/2);
                               
                               innercount =innercount+1;
                               seeddag = tempdag;
                               allscore(j)= temp_score;
                               SampleDistribution{j} = temp_gamma;
                               hiddenGraph_Ps{j} = temp_hiddenGraph_Ps;
                           
                           end;

                   end;
               end;


            end; % end for j
        end; % end for i

        if innercount == 0
            outcount = outcount +1;
        end;

    end;  % end while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % end function