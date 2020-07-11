function output = HPS2(filename,type,pow,win_size,hop_size,medfilt_len,outfile)
%harmonic/percussive separation processing
%function output1 = HPS2(filename,type,pow,win_size,hop_size,medfilt_len,outfile)

%read the signal and make mono
[wave,fs]=audioread(filename);
wave = wave(:,1);

%sine window
win = sin(pi/win_size*(1:win_size)');

%buffer, take fft, and magnitude
y = buffer(wave,win_size,win_size-hop_size);
y = bsxfun(@times,win,y);
Y = fft(y,[],1);
[th,r] = cart2pol(real(Y),imag(Y));
PS = zeros(size(r,1)/2,size(r,2));
HS = zeros(size(r,1)/2,size(r,2));
out = zeros(size(r,1)/2,size(r,2));

%median filter lengths
if mod(medfilt_len,2) == 1
    
    L1 = ((medfilt_len-1)/2);
    L2 = ((medfilt_len-1)/2);
    
elseif mod(medfilt_len,2) == 0
   
    L1 = ceil((medfilt_len-1)/2);
    L2 = floor((medfilt_len-1)/2);
    
end

%median filter percussive suppression mask
for n = 1:size(r,1)/2

    for m = 1:L1
    
        PS(n,m) = median(r(n,1:m+L2));

    end
    
    for m = L1+1:size(r,2)-L2

        PS(n,m) = median(r(n,m-L1:m+L2));

    end

	for m = size(r,2)-L2+1:size(r,2)

        PS(n,m) = median(r(n,m-L1:size(r,2)));

    end

end

%median filter harmonic suppression mask
for m = 1:size(r,2)

    for n = 1:L1

        HS(n,m) = median(r(1:n+L2,m));

    end

    for n = L1+1:size(r,1)/2-L2

        HS(n,m) = median(r(n-L1:n+L2,m));

    end

    for n = size(r,1)/2-L2+1:size(r,1)/2

        HS(n,m) = median(r(n-L1:size(r,1)/2,m));

    end
    
end

%apply Weiner mask
switch type
    
    %percussive mask
    case 'harm'
        
        for n = 1:size(r,1)/2
            
            for m = 1:size(r,2)
                
                out(n,m) = power(PS(n,m),pow)/(power(HS(n,m),pow) + power(PS(n,m),pow)) * r(n,m);
                
            end
            
        end
    
    %harmonic mask
    case 'perc'
        
        for n = 1:size(r,1)/2
            
            for m = 1:size(r,2)
                
                out(n,m) = power(HS(n,m),pow)/(power(HS(n,m),pow)+power(PS(n,m),pow)) * r(n,m);
            
            end
            
        end
        
end

%mirror until length of fft
R = [out(1:end,:);out(end:-1:1,:)];
TH = [th(1:end/2,:);th(end/2:-1:1,:)];

%polar to cartesian
[re, im] = pol2cart(TH,R);
output = real(ifft(re + i*im,[],1));

%window
output = bsxfun(@times,win,output);

%allocate output
output1 = zeros(size(output,2)*size(output,1)/2-1,2);

%overlap add
starting_samp = 1; 
for m = 1:size(output,2)-1
    
    output1(starting_samp:starting_samp + win_size-1) = output1(starting_samp:starting_samp + win_size-1) + output(:,m)';
    
    starting_samp = starting_samp + hop_size;

end

%output
output1(:,2) = output1(:,1);

%output the result
audiowrite(outfile,output1,fs);

end