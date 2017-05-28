

function AddHeaders(savename, titlenames, MatrixtoGet)
% AddHeaders save given matrix with the given header names, as a .csv file
% AddHeaders(savename, titlenames, MatrixtoGet)
%{  savename is string
%   titlenames is a cell array of headers for the output
%   MatrixtoGet is the array to be headered and saved
%}


% fileID = fopen(savename,'w+');
%     fprintf(fileID, headerwhere, titles);
%     fprintf(fileID, datawhere, MatrixtoGet');
% fclose(fileID);




% These are the formatSpec parameters 
headerpostions = {'%1s','%2s','%3s','%4s','%5s','%6s','%7s','%8s','%9s',...
    '%10s','%11s','%12s','%13s','%14s','%15s','%16s','%17s','%18s','%19s',...
    '%20s','%21s','%22s','%23s','%24s','%25s','%26s','%27s','%28s','%29s',...
    '%30s','%31s','%32s','%33s','%34s','%35s','%36s','%37s','%38s','%39s',...
    '%40s','%41s','%42s','%43s','%44s','%45s','%46s','%47s','%48s','%49s',...
    '%50s','%51s','%52s','%53s','%54s','%55s','%56s','%57s','%58s','%59s',...
    '%60s','%61s','%62s','%63s','%64s','%65s','%66s','%67s','%68s','%69s',...
    '%70s','%71s','%72s','%73s','%74s','%75s','%76s','%77s','%78s','%79s',...
    '%80s','%81s','%82s','%83s','%84s','%85s','%86s','%87s','%88s','%89s',...
    '%90s','%91s','%92s','%93s','%94s','%95s','%96s','%97s','%98s','%99s',...
    '%100s','%101s','%102s','%103s','%104s','%105s','%106s','%107s','%108s','%109s',...
    '%110s','%111s','%112s','%113s','%114s','%115s','%116s','%117s','%118s','%119s',...
    '%120s','%121s','%122s','%123s','%124s','%125s','%126s','%127s','%128s','%129s',...
    '%130s','%131s','%132s','%133s','%134s','%135s','%136s','%137s','%138s','%139s',...
    '%140s','%141s','%142s','%143s','%144s','%145s','%146s','%147s','%148s','%149s',...
    '%150s','%151s','%152s','%153s','%154s','%155s','%156s','%157s','%158s','%159s',...
    '%160s','%161s','%162s','%163s','%164s','%165s','%166s','%167s','%168s','%169s',...
    '%170s','%171s','%172s','%173s','%174s','%175s','%176s','%177s','%178s','%179s',...
    '%180s','%181s','%182s','%183s','%184s','%185s','%186s','%187s','%188s','%189s',...
    '%190s','%191s','%192s','%193s','%194s','%195s','%196s','%197s','%198s','%199s',...
    '%200s','%201s','%202s','%203s','%204s','%205s','%206s','%207s','%208s','%209s',...
    '%210s','%211s','%212s','%213s','%214s','%215s','%216s','%217s','%218s','%219s',...
    '%220s','%221s','%222s','%223s','%224s','%225s','%226s','%227s','%228s','%229s',...
    '%230s','%231s','%232s','%233s','%234s','%235s','%236s','%237s','%238s','%239s',...
    '%240s','%241s','%242s','%243s','%244s','%245s','%246s','%247s','%248s','%249s',...
    '%250s','%251s','%252s','%253s','%254s','%255s','%256s','%257s','%258s','%259s',...
    '%260s','%261s','%262s','%263s','%264s','%265s','%266s','%267s','%268s','%269s',...
    '%270s','%271s','%272s','%273s','%274s','%275s','%276s','%277s','%278s','%279s',...
    '%280s','%281s','%282s','%283s','%284s','%285s','%286s','%287s','%288s','%289s',...
    '%290s','%291s','%292s','%293s','%294s','%295s','%296s','%297s','%298s','%299s',...
    '%300s','%301s','%302s','%303s','%304s','%305s','%306s','%307s','%308s','%309s',...
    '%310s','%311s','%312s','%313s','%314s','%315s','%316s','%317s','%318s','%319s',...
    '%320s','%321s','%322s','%323s','%324s','%325s','%326s','%327s','%328s','%329s',...
    '%330s','%331s','%332s','%333s','%334s','%335s','%336s','%337s','%338s','%339s',...
    '%340s','%341s','%342s','%343s','%344s','%345s','%346s','%347s','%348s','%349s',...
    '%350s','%351s','%352s','%353s','%354s','%355s','%356s','%357s','%358s','%359s',...
    '%360s','%361s','%362s','%363s','%364s','%365s','%366s','%367s','%368s','%369s',...
    '%370s','%371s','%372s','%373s','%374s','%375s','%376s','%377s','%378s','%379s',...
    '%380s','%381s','%382s','%383s','%384s','%385s','%386s','%387s','%388s','%389s',...
    '%390s','%391s','%392s','%393s','%394s','%395s','%396s','%397s','%398s','%399s',...
    '%400s','%401s','%402s','%403s','%404s','%405s','%406s','%407s','%408s','%409s',...
    '%410s','%411s','%412s','%413s','%414s','%415s','%416s','%417s','%418s','%419s',...
    '%420s','%421s','%422s','%423s','%424s','%425s','%426s','%427s','%428s','%429s',...
    '%430s','%431s','%432s','%433s','%434s','%435s','%436s','%437s','%438s','%439s',...
    '%440s','%441s','%442s','%443s','%444s','%445s','%446s','%447s','%448s','%449s',...
    '%450s','%451s','%452s','%453s','%454s','%455s','%456s','%457s','%458s','%459s',...
    '%460s','%461s','%462s','%463s','%464s','%465s','%466s','%467s','%468s','%469s',...
    '%470s','%471s','%472s','%473s','%474s','%475s','%476s','%477s','%478s','%479s',...
    '%480s','%481s','%482s','%483s','%484s','%485s','%486s','%487s','%488s','%489s',...
    '%490s','%491s','%492s','%493s','%494s','%495s','%496s','%497s','%498s','%499s'};

datapostions = {'%1g','%2g','%3g','%4g','%5g','%6g','%7g','%8g','%9g',...
    '%10g','%11g','%12g','%13g','%14g','%15g','%16g','%17g','%18g','%19g',...
    '%20g','%21g','%22g','%23g','%24g','%25g','%26g','%27g','%28g','%29g',...
    '%30g','%31g','%32g','%33g','%34g','%35g','%36g','%37g','%38g','%39g',...
    '%40g','%41g','%42g','%43g','%44g','%45g','%46g','%47g','%48g','%49g',...
    '%50g','%51g','%52g','%53g','%54g','%55g','%56g','%57g','%58g','%59g',...
    '%60g','%61g','%62g','%63g','%64g','%65g','%66g','%67g','%68g','%69g',...
    '%70g','%71g','%72g','%73g','%74g','%75g','%76g','%77g','%78g','%79g',...
    '%80g','%81g','%82g','%83g','%84g','%85g','%86g','%87g','%88g','%89g',...
    '%90s','%91s','%92s','%93s','%94s','%95s','%96s','%97s','%98s','%99s',...
    '%100g','%101g','%102g','%103g','%104g','%105g','%106g','%107g','%108g','%109g',...
    '%110g','%111g','%112g','%113g','%114g','%115g','%116g','%117g','%118g','%119g',...
    '%120g','%121g','%122g','%123g','%124g','%125g','%126g','%127g','%128g','%129g',...
    '%130g','%131g','%132g','%133g','%134g','%135g','%136g','%137g','%138g','%139g',...
    '%140g','%141g','%142g','%143g','%144g','%145g','%146g','%147g','%148g','%149g',...
    '%150g','%151g','%152g','%153g','%154g','%155g','%156g','%157g','%158g','%159g',...
    '%160g','%161g','%162g','%163g','%164g','%165g','%166g','%167g','%168g','%169g',...
    '%170g','%171g','%172g','%173g','%174g','%175g','%176g','%177g','%178g','%179g',...
    '%180g','%181g','%182g','%183g','%184g','%185g','%186g','%187g','%188g','%189g',...
    '%190g','%191g','%192g','%193g','%194g','%195g','%196g','%197g','%198g','%199g',...
    '%200g','%201g','%202g','%203g','%204g','%205g','%206g','%207g','%208g','%209g',...
    '%210g','%211g','%212g','%213g','%214g','%215g','%216g','%217g','%218g','%219g',...
    '%220g','%221g','%222g','%223g','%224g','%225g','%226g','%227g','%228g','%229g',...
    '%230g','%231g','%232g','%233g','%234g','%235g','%236g','%237g','%238g','%239g',...
    '%240g','%241g','%242g','%243g','%244g','%245g','%246g','%247g','%248g','%249g',...
    '%250g','%251g','%252g','%253g','%254g','%255g','%256g','%257g','%258g','%259g',...
    '%260g','%261g','%262g','%263g','%264g','%265g','%266g','%267g','%268g','%269g',...
    '%270g','%271g','%272g','%273g','%274g','%275g','%276g','%277g','%278g','%279g',...
    '%280g','%281g','%282g','%283g','%284g','%285g','%286g','%287g','%288g','%289g',...
    '%290g','%291g','%292g','%293g','%294g','%295g','%296g','%297g','%298g','%299g',...
    '%300g','%301g','%302g','%303g','%304g','%305g','%306g','%307g','%308g','%309g',...
    '%310g','%311g','%312g','%313g','%314g','%315g','%316g','%317g','%318g','%319g',...
    '%320g','%321g','%322g','%323g','%324g','%325g','%326g','%327g','%328g','%329g',...
    '%330g','%331g','%332g','%333g','%334g','%335g','%336g','%337g','%338g','%339g',...
    '%340g','%341g','%342g','%343g','%344g','%345g','%346g','%347g','%348g','%349g',...
    '%350g','%351g','%352g','%353g','%354g','%355g','%356g','%357g','%358g','%359g',...
    '%360g','%361g','%362g','%363g','%364g','%365g','%366g','%367g','%368g','%369g',...
    '%370g','%371g','%372g','%373g','%374g','%375g','%376g','%377g','%378g','%379g',...
    '%380g','%381g','%382g','%383g','%384g','%385g','%386g','%387g','%388g','%389g',...
    '%390g','%391g','%392g','%393g','%394g','%395g','%396g','%397g','%398g','%399g',...
    '%400g','%401g','%402g','%403g','%404g','%405g','%406g','%407g','%408g','%409g',...
    '%410g','%411g','%412g','%413g','%414g','%415g','%416g','%417g','%418g','%419g',...
    '%420g','%421g','%422g','%423g','%424g','%425g','%426g','%427g','%428g','%429g',...
    '%430g','%431g','%432g','%433g','%434g','%435g','%436g','%437g','%438g','%439g',...
    '%440g','%441g','%442g','%443g','%444g','%445g','%446g','%447g','%448g','%449g',...
    '%450g','%451g','%452g','%453g','%454g','%455g','%456g','%457g','%458g','%459g',...
    '%460g','%461g','%462g','%463g','%464g','%465g','%466g','%467g','%468g','%469g',...
    '%470g','%471g','%472g','%473g','%474g','%475g','%476g','%477g','%478g','%479g',...
    '%480g','%481g','%482g','%483g','%484g','%485g','%486g','%487g','%488g','%489g',...
    '%490g','%491g','%492g','%493g','%494g','%495g','%496g','%497g','%498g','%499g'};


% Get only the number of placer strings as there are data columns in the matrix

% Make the cell array of title strings so they can be written in
titles = strjoin(titlenames(1:size(MatrixtoGet,2)),',');

fileID = fopen(savename,'w+');

    headerFormat = strcat(strjoin(headerpostions(1:(size(MatrixtoGet,2))),','),'\n');
    fprintf(fileID, headerFormat, titles);

    dataFormat = strcat('\n',strjoin(datapostions(1:(size(MatrixtoGet,2))),','));
    for row = 1:size(MatrixtoGet,1),
        fprintf(fileID, dataFormat, MatrixtoGet(row,:));
    end
    
fclose(fileID);
% For some reason, the default is to write the matrix by rows, not by colunmns... hence, the ' inversion operator

end


