"""Capture only this application's Linux X11 window for visual verification."""
import ctypes
import argparse
import json
import os
import struct
import sys
import time
import tkinter as tk
import zlib
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from app import ProteinLocApp
from pipeline import list_jobs,read_json


class XImage(ctypes.Structure):
    _fields_=[('width',ctypes.c_int),('height',ctypes.c_int),('xoffset',ctypes.c_int),('format',ctypes.c_int),
              ('data',ctypes.c_void_p),('byte_order',ctypes.c_int),('bitmap_unit',ctypes.c_int),
              ('bitmap_bit_order',ctypes.c_int),('bitmap_pad',ctypes.c_int),('depth',ctypes.c_int),
              ('bytes_per_line',ctypes.c_int),('bits_per_pixel',ctypes.c_int),
              ('red_mask',ctypes.c_ulong),('green_mask',ctypes.c_ulong),('blue_mask',ctypes.c_ulong),
              ('obdata',ctypes.c_void_p),('functions',ctypes.c_void_p*6)]


def capture(root,path):
    x=ctypes.CDLL('libX11.so.6')
    x.XOpenDisplay.argtypes=[ctypes.c_char_p]
    x.XOpenDisplay.restype=ctypes.c_void_p
    x.XGetImage.argtypes=[ctypes.c_void_p,ctypes.c_ulong,ctypes.c_int,ctypes.c_int,ctypes.c_uint,ctypes.c_uint,ctypes.c_ulong,ctypes.c_int]
    x.XGetImage.restype=ctypes.POINTER(XImage)
    x.XDestroyImage.argtypes=[ctypes.POINTER(XImage)]
    x.XCloseDisplay.argtypes=[ctypes.c_void_p]
    display=x.XOpenDisplay(os.environ['DISPLAY'].encode())
    if not display:
        raise RuntimeError('Cannot open Linux display')
    image=x.XGetImage(display,root.winfo_id(),0,0,root.winfo_width(),root.winfo_height(),ctypes.c_ulong(-1).value,2)
    if not image:
        raise RuntimeError('Cannot read application window image')
    try:
        info=image.contents
        if info.bits_per_pixel != 32 or info.byte_order != 0 or info.red_mask != 0xff0000:
            raise RuntimeError('Unexpected X11 pixel format')
        raw=ctypes.string_at(info.data,info.bytes_per_line*info.height)
        bgra=np.frombuffer(raw,dtype=np.uint8).reshape(info.height,info.bytes_per_line)[:,:info.width*4].reshape(info.height,info.width,4)
        rgb=bgra[:,:,[2,1,0]]
        data=b''.join(b'\x00'+row.tobytes() for row in rgb)
        def chunk(kind,value):
            return struct.pack('>I',len(value))+kind+value+struct.pack('>I',zlib.crc32(kind+value)&0xffffffff)
        path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',info.width,info.height,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(data,5))+chunk(b'IEND',b''))
    finally:
        x.XDestroyImage(image)
        x.XCloseDisplay(display)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--job',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/gui.png')
    args=parser.parse_args()
    root=tk.Tk()
    app=ProteinLocApp(root,runs_dir=ROOT/'validation/runs')
    root.geometry('1300x880')
    references=[s for s in list_jobs(ROOT/'validation/runs') if s['count']==3 and s['status'] in ('completed','partial')]
    if args.job:
        app.load_job(args.job.resolve())
        group=read_json(args.job/'job.json')['organism']
        from app import GROUPS
        app.group.set(next(label for label,value in GROUPS.items() if value==group))
    elif references:
        app.load_job(references[0]['path'])
    started=time.monotonic()
    while time.monotonic()-started<5:
        root.update()
        time.sleep(.05)
    children=app.tree.get_children()
    if children:
        app.tree.selection_set(children[0])
        app.show_detail()
    root.update()
    capture(root,args.output)
    print(json.dumps({'window':root.title(),'geometry':root.geometry(),'shown_records':len(app.records),'screenshot':str(args.output)},ensure_ascii=False))
    root.destroy()


if __name__=='__main__':
    main()
