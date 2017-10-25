function cost = corrcheck_PBC(L,Jmaxpos_vec,Si,Sj,chi_inc)
%cost = corrcheck_PBC(L,Jmaxpos_vec,Si,Sj,chi_inc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dMERA - corrcheck_PBC
% finds cost of contraction that would be done using corr_PBC
% 
% Andrew Goldsborough - 16/01/2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keep track of cost
cost = 0;

%flag for two sites having combined and size of operator
if min(PBC_pos(Si-Sj,L),PBC_pos(Sj-Si,L)) == 1
    %neighbours are trivially successful
    combined = 2;
else
    combined = 0;
end

%check left and right are correct
if Si > Sj 
    %swap L and R
    temp = Sj;
    Sj = Si;
    Si = temp;
end

%run over algorithm
for iteration = 1:(L/2)-2
    
    %store current length and position
    Lcur = L-2*(iteration-1);
    Jmaxpos = Jmaxpos_vec(iteration);
    
    %check combined
    if combined == 0
        %%% not combined %%%
        
        %%% update both together %%%
        if min(mod(Si-Jmaxpos,Lcur),mod(Jmaxpos-Si,Lcur)) <= 2 && min(mod(Sj-Jmaxpos,Lcur),mod(Jmaxpos-Sj,Lcur)) <= 2
            
            %check left and right are correct
            if PBC_pos(Sj-Jmaxpos+3,Lcur) < PBC_pos(Si-Jmaxpos+3,Lcur)
                %swap L and R
                temp = Sj;
                Sj = Si;
                Si = temp;
            end
            
            if Si == PBC_pos(Jmaxpos-1,Lcur) && Sj == PBC_pos(Jmaxpos+1,Lcur)
                %store as L
                if chi_inc(iteration) == 0
                    %corr cost chi^6
                    cost = max(cost,6);
                else %chi_inc(iteration) == 1
                    %corr_inc cost chi^6
                    cost = max(cost,6);
                end
                
                %remove two sites
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 1,Lcur-2);
                end
                
                Sj = 0;
                combined = 2;
                
            elseif Si == PBC_pos(Jmaxpos-2,Lcur) && Sj == PBC_pos(Jmaxpos+2,Lcur)
                if chi_inc(iteration) == 0
                    %4-site operator
                    
                    %make4 cost chi^10
                    cost = max(cost,10);
                    
                    %remove two sites set 4-flag
                    if Jmaxpos == Lcur
                        Si = PBC_pos(Jmaxpos - 3,Lcur-2);
                    else
                        Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                    end
                    
                    Sj = 0;
                    combined = 4;
                    
                else %chi_inc(iteration) == 1
                    %splits into two
                    
                    %h_1 and h_5 cost chi^6
                    cost = max(cost,6);
                    
                    %remove two sites
                    if Jmaxpos == Lcur
                        Si = PBC_pos(Jmaxpos - 3,Lcur-2);
                        Sj = 1;
                    else
                        Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                        Sj = PBC_pos(Jmaxpos,Lcur-2);
                    end
                end
                
            elseif Si == PBC_pos(Jmaxpos-2,Lcur) && Sj == PBC_pos(Jmaxpos,Lcur)
                %creates a 3-site operator (1)
                if chi_inc(iteration) == 0
                    
                    %make3_1 cost chi^8
                    cost = max(cost,8);
                    
                else %chi_inc(iteration) == 1
                    %make3_1_inc cost chi^8
                    cost = max(cost,8);
                end
                
                %remove two sites, set 3 flag
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 3,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                end
                
                Sj = 0;
                combined = 3;
                
            elseif Si == PBC_pos(Jmaxpos-2,Lcur) && Sj == PBC_pos(Jmaxpos+1,Lcur)
                %creates a 3-site operator (2)
                if chi_inc(iteration) == 0
                    
                    %make3_2 cost chi^8
                    cost = max(cost,8);
                    
                else %chi_inc(iteration) == 1
                    %make3_2_inc cost chi^8
                    cost = max(cost,8);
                end
                
                %remove two sites, set 3 flag
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 3,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                end
                
                Sj = 0;
                combined = 3;
                
            elseif Si == PBC_pos(Jmaxpos-1,Lcur) && Sj == PBC_pos(Jmaxpos+2,Lcur)
                %creates a 3-site operator (3)
                if chi_inc(iteration) == 0
                    %make3_3 cost chi^8
                    cost = max(cost,8);
                    
                else %chi_inc(iteration) == 1
                    %make3_3 cost chi^8
                    cost = max(cost,8);
                end
                
                %remove two sites, set 3 flag
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 1,Lcur-2);
                end
                
                Sj = 0;
                combined = 3;
                
            elseif Si == PBC_pos(Jmaxpos,Lcur) && Sj == PBC_pos(Jmaxpos+2,Lcur)
                %creates a 3-site operator (4)
                if chi_inc(iteration) == 0
                    %make3_4 cost chi^8
                    cost = max(cost,8);
                    
                else %chi_inc(iteration) == 1
                    %make3_4_inc cost chi^8
                    cost = max(cost,8);
                end
                
                %remove two sites, set 3 flag
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 1,Lcur-2);
                end
                
                Sj = 0;
                combined = 3;
            end
            
            %%% just update SzL %%%
        elseif Si == PBC_pos(Jmaxpos - 2,Lcur)
            %diagram 1: leftmost
            if chi_inc(iteration) == 0
                %h_1 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_1_inc cost chi^6
                cost = max(cost,6);
            end
            
            Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            Sj = PBC_pos(Sj - 2,Lcur-2);
        elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
            %diagram 2: left
            if chi_inc(iteration) == 0
                %h_2 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_2_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            Sj = PBC_pos(Sj - 2,Lcur-2);
        elseif Si == Jmaxpos
            %diagram 3: centre
            if chi_inc(iteration) == 0
                %h_3 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_3_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            Sj = PBC_pos(Sj - 2,Lcur-2);
        elseif Si == PBC_pos(Jmaxpos + 1,Lcur)
            %diagram 4: right
            if chi_inc(iteration) == 0
                %h_4 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_4_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos-2,Lcur-2);
                Sj = PBC_pos(Sj - 1,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos-1,Lcur-2);
                Sj = PBC_pos(Sj - 2,Lcur-2);
            end
            
        elseif Si == PBC_pos(Jmaxpos + 2,Lcur)
            %diagram 5 rightmost
            if chi_inc(iteration) == 0
                %h_5 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_5_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos-1,Lcur-2);
                Sj = PBC_pos(Sj - 1,Lcur-2);
            elseif Jmaxpos == Lcur-1
                %no change in indices
            else
                Si = PBC_pos(Jmaxpos,Lcur-2);
                Sj = PBC_pos(Sj - 2,Lcur-2);
            end
            
            %%% just update SzR %%%
        elseif Sj == PBC_pos(Jmaxpos - 2,Lcur)
            %diagram 1: leftmost
            if chi_inc(iteration) == 0
                %h_1 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_1_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Si - 1,Lcur-2);
                Sj = PBC_pos(Jmaxpos - 3,Lcur-2);
            elseif Jmaxpos == 1
                Si = PBC_pos(Si - 2,Lcur-2);
                Sj = Lcur - 3;
            elseif Jmaxpos == 2
                Si = PBC_pos(Si - 2,Lcur-2);
                Sj = Lcur - 2;
            else
                Sj = PBC_pos(Jmaxpos - 2,Lcur-2);
            end
            
        elseif Sj == PBC_pos(Jmaxpos - 1,Lcur)
            %diagram 2: left
            if chi_inc(iteration) == 0
                %h_2 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_2_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Si - 1,Lcur-2);
                Sj = PBC_pos(Jmaxpos - 2,Lcur-2);
            elseif Jmaxpos == 1
                Si = PBC_pos(Si - 2,Lcur-2);
                Sj = Lcur - 2;
            else
                Sj = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
        elseif Sj == Jmaxpos
            %diagram 3: centre
            if chi_inc(iteration) == 0
                %h_3 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_3_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Si - 1,Lcur-2);
                Sj = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Sj = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
        elseif Sj == PBC_pos(Jmaxpos + 1,Lcur)
            %diagram 4: right
            if chi_inc(iteration) == 0
                %h_4 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_4_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            Sj = PBC_pos(Jmaxpos - 1,Lcur-2);
        elseif Sj == PBC_pos(Jmaxpos + 2,Lcur)
            %diagram 5 rightmost
            if chi_inc(iteration) == 0
                %h_5 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_5_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            Sj = PBC_pos(Jmaxpos,Lcur-2);
        else
            %no update
            if Jmaxpos < Si
                %two sites removed from the left of both sites
                Si = Si - 2;
                Sj = Sj - 2;
            elseif Jmaxpos == Lcur
                %remove one site from the left one from the right
                Si = Si - 1;
                Sj = Sj - 1;
            elseif Jmaxpos > Sj;
                %to the right of both sites => no action
            else
                Sj = Sj - 2;
            end
        end
        
        %check if L and R have swapped
        if Si > Sj && Sj ~= 0
            temp = Sj;
            Sj = Si;
            Si = temp;
        end
        
    elseif combined == 2
        %%% combined, 2 site %%%
        
        %just update SzL
        if Si == PBC_pos(Jmaxpos - 2,Lcur)
            %diagram 1: leftmost
            if chi_inc(iteration) == 0
                %h_1 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_1_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 3,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
            %diagram 2: left
            if chi_inc(iteration) == 0
                %h_2 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_2_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
        elseif Si == Jmaxpos
            %diagram 3: centre
            if chi_inc(iteration) == 0
                %h_3 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_3_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos + 1,Lcur)
            %diagram 4: right
            if chi_inc(iteration) == 0
                %h_4 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %h_4_inc cost chi^7
                cost = max(cost,7);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos + 2,Lcur)
            %diagram 5 rightmost
            if chi_inc(iteration) == 0
                %h_5 cost chi^6
                cost = max(cost,6);
            else %chi_inc(iteration) == 1
                %h_5_inc cost chi^6
                cost = max(cost,6);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos,Lcur-2);
            end
        else
            %no update
            if Jmaxpos < Si
                %two sites removed from the left of both sites
                Si = Si - 2;
            elseif Jmaxpos == Lcur
                %remove one site from the left one from the right
                Si = Si - 1;
                
                %to the right of both sites => no action
            end
        end
    elseif combined == 3
        %%% combined, 3 site %%%
        
        %just update SzL
        if Si == PBC_pos(Jmaxpos - 3,Lcur)
            %propagate 3 (1)
            if chi_inc(iteration) == 0
                %prop3_1 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %prop3_1_inc cost chi^8
                cost = max(cost,8);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 4,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 3,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos - 2,Lcur)
            %propagate 3 (2)
            if chi_inc(iteration) == 0
                %prop3_2 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop3_2_inc cost chi^9
                cost = max(cost,9);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 3,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
            %propagate 3 (3)
            if chi_inc(iteration) == 0
                %prop3_3 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %prop3_3_inc cost chi^8
                cost = max(cost,8);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
            %now a two-site operator
            combined = 2;
        elseif Si == Jmaxpos
            %propagate 3 (4)
            if chi_inc(iteration) == 0
                %prop3_4 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %prop3_4_inc cost chi^8
                cost = max(cost,8);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
            %now a two-site operator
            combined = 2;
        elseif Si == PBC_pos(Jmaxpos + 1,Lcur)
            %propagate 3 (5)
            if chi_inc(iteration) == 0
                %prop3_5 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop3_5_inc cost chi^9
                cost = max(cost,9);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
        elseif Si == PBC_pos(Jmaxpos + 2,Lcur)
            %propagate 3 (6)
            if chi_inc(iteration) == 0
                %prop3_6 cost chi^8
                cost = max(cost,8);
            else %chi_inc(iteration) == 1
                %prop3_6_inc cost chi^8
                cost = max(cost,8);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos,Lcur-2);
            end
        else
            %no update
            if Jmaxpos < Si
                %two sites removed from the left of both sites
                Si = Si - 2;
            elseif Jmaxpos == Lcur
                %remove one site from the left one from the right
                Si = Si - 1;
                
                %to the right of both sites => no action
            end
        end
    else %combined == 4
        %%% combined, 4 site %%%
        
        %just update SzL
        if Si == PBC_pos(Jmaxpos - 4,Lcur)
            if Lcur == 6
                %PBC term (6)
                if chi_inc(iteration) == 0
                    %prop4_PBC6 cost chi^12
                    cost = max(cost,12);
                else %chi_inc(iteration) == 1
                    %prop4_PBC6_inc cost chi^10
                    cost = max(cost,10);
                end
                
                %remove two sites
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 3,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 2,Lcur-2);
                end
            else
                %propagate 4 (1)
                if chi_inc(iteration) == 0
                    %prop4_1 cost chi^10
                    cost = max(cost,10);
                else %chi_inc(iteration) == 1
                    %prop4_1_inc cost chi^10
                    cost = max(cost,10);
                end
                
                %remove two sites
                if Jmaxpos == Lcur
                    Si = PBC_pos(Jmaxpos - 5,Lcur-2);
                else
                    Si = PBC_pos(Jmaxpos - 4,Lcur-2);
                end
            end
            
        elseif Si == PBC_pos(Jmaxpos - 3,Lcur)
            %propagate 4 (2)
            if chi_inc(iteration) == 0
                %prop4_2 cost chi^12
                cost = max(cost,12);
            else %chi_inc(iteration) == 1
                %prop4_2_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 4,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 3,Lcur-2);
            end
            
        elseif Si == PBC_pos(Jmaxpos - 2,Lcur)
            %propagate 4 (3)
            if chi_inc(iteration) == 0
                %prop4_3 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop4_3_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 3,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            end
            
            %now a three-site operator
            combined = 3;
            
        elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
            %propagate 4 (4)
            if chi_inc(iteration) == 0
                %prop4_4 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop4_4_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
            %now a two-site operator
            combined = 2;
            
        elseif Si == PBC_pos(Jmaxpos,Lcur)
            %propagate 4 (5)
            if chi_inc(iteration) == 0
                %prop4_5 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop4_5_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
            %now a three-site operator
            combined = 3;
            
        elseif Si == PBC_pos(Jmaxpos + 1,Lcur)
            %propagate 4 (6)
            if chi_inc(iteration) == 0
                %prop4_6 cost chi^12
                cost = max(cost,12);
            else %chi_inc(iteration) == 1
                %prop4_6_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 2,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            end
            
        elseif Si == PBC_pos(Jmaxpos + 2,Lcur)
            %propagate 4 (7)
            if chi_inc(iteration) == 0
                %prop4_7 cost chi^10
                cost = max(cost,10);
            else %chi_inc(iteration) == 1
                %prop4_7_inc cost chi^10
                cost = max(cost,10);
            end
            
            %remove two sites
            if Jmaxpos == Lcur
                Si = PBC_pos(Jmaxpos - 1,Lcur-2);
            else
                Si = PBC_pos(Jmaxpos,Lcur-2);
            end
        else
            %no update
            if Jmaxpos < Si
                %two sites removed from the left of both sites
                Si = Si - 2;
            elseif Jmaxpos == Lcur
                %remove one site from the left one from the right
                Si = Si - 1;
                
                %to the right of both sites => no action
            end
        end
    end
end

%alternative top
iteration = (L/2)-1;
Jmaxpos = Jmaxpos_vec(iteration);
Lcur = 4;

%check combined
if combined == 0
    %%% not combined %%%
    
    %check left and right are correct
    if PBC_pos(Sj-Jmaxpos+2,Lcur) < PBC_pos(Si-Jmaxpos+2,Lcur)
        %swap L and R
        temp = Sj;
        Sj = Si;
        Si = temp;
    end
    
    if Si == PBC_pos(Jmaxpos-1,Lcur) && Sj == PBC_pos(Jmaxpos+1,Lcur)
        %corr_PBC1_at cost chi^6
        cost = max(cost,6);
    elseif Si == PBC_pos(Jmaxpos,Lcur) && Sj == PBC_pos(Jmaxpos+2,Lcur)    
        %corr_PBC2_at cost chi^6
        cost = max(cost,6);
    end

elseif combined == 2
    %%% combined, 2 site %%%
    if Si == PBC_pos(Jmaxpos - 2,Lcur)
        %PBC term
        %h_PBC1_at cost chi^6
        cost = max(cost,6);
    elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
        %diagram 2: left
        %h_PBC2_at cost chi^6
        cost = max(cost,6);
    elseif Si == Jmaxpos
        %diagram 3: centre
        %h_PBC3_at cost chi^6
        cost = max(cost,6);        
    elseif Si == PBC_pos(Jmaxpos + 1,Lcur)
        %diagram 4: right
        %h_PBC4_at cost chi^6
        cost = max(cost,6);
    end
    
elseif combined == 3
    %%% combined, 3 site %%%
    if Si == PBC_pos(Jmaxpos - 3,Lcur)
        %PBC term (1)
        %prop3_PBC1_at cost chi^7
        cost = max(cost,7);
    elseif Si == PBC_pos(Jmaxpos - 2,Lcur)
        %PBC term (2)
        %prop3_PBC2_at cost chi^7
        cost = max(cost,7);
    elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
        %propagate 3 (3)
        %prop3_PBC3_at cost chi^7
        cost = max(cost,7);
    elseif Si == Jmaxpos
        %propagate 3 (4)
        %prop3_PBC4_at cost chi^7
        cost = max(cost,7);
    end
    
else %combined == 4
    %%% combined, 4 site %%%
    if Si == PBC_pos(Jmaxpos,Lcur)
        %PBC term (1)
        %prop4_PBC1_at cost chi^8
        cost = max(cost,8);
    elseif Si == PBC_pos(Jmaxpos - 3,Lcur)
        %PBC term (2)
        %prop4_PBC2_at cost chi^8
        cost = max(cost,8);
    elseif Si == PBC_pos(Jmaxpos - 2,Lcur)
        %PBC term (3)
        %prop4_PBC3_at cost chi^8
        cost = max(cost,8);
    elseif Si == PBC_pos(Jmaxpos - 1,Lcur)
        %propagate 4 (4)
        %prop4_PBC4_at cost chi^8
        cost = max(cost,8);
    end
end
end