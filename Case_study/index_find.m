function[ls] = index_find(ls)
        if ls <= 10
            combo = nchoosek(1:10, 1);
            ls = combo(ls, :);
        end
        if ls > 10 & ls <= 55
            combo = nchoosek(1:10, 2);
            ls = combo(ls - 10, :);
        end
        if ls > 55 & ls <= 175
            combo = nchoosek(1:10, 3);
            ls = combo(ls - 55, :);
        end
        if ls > 175 & ls <= 385
            combo = nchoosek(1:10, 4);
            ls = combo(ls - 175, :);
        end
        if ls > 385 & ls <= 637
            combo = nchoosek(1:10, 5);
            ls = combo(ls - 385, :);
        end
        if ls > 637 & ls <= 847
            combo = nchoosek(1:10, 6);
            ls = combo(ls - 637, :);
        end
        if ls > 847 & ls <= 967
            combo = nchoosek(1:10, 7);
            ls = combo(ls - 847, :);
        end
        if ls > 967 & ls <= 1012
            combo = nchoosek(1:10, 8);
            ls = combo(ls - 967, :);
        end
        if ls > 1012 & ls <= 1022
            combo = nchoosek(1:10, 9);
            ls = combo(ls - 1012, :);
        end
        if ls == 1023
            combo = nchoosek(1:10, 10);
            ls = combo(ls - 1022, :);
        end
end