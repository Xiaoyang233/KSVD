function [Pic,PathName,PicName] = loadImage()

   [PicName,PathName,~] = uigetfile('d:\m_sou_file\����ʦ��ң������\*.*','Selcet a file');
   if(PicName == 0)
        return;
    end
    PicStr = [PathName,PicName];
    Pic = imread(PicStr);
    if(length(size(Pic))==3)
        Pic = rgb2gray(Pic);
    end
   
end