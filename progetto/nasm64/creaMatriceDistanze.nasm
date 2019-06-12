%include "sseutils64.nasm"

extern distanza

section .data
d1 db '%ld ', 10, 0
d3 db '%ld %ld %ld ', 10, 0
d2 db '%ld %ld ', 10, 0
;f db '%f ', 10, 0
break db 'breakpoint', 10, 0
section .bss
section .text

global creaMatriceDistanze

; codebook equ 12   rsi
; input equ 8       rdi

d equ 28
m equ 40
k equ 44
distanze_simmetriche equ 128

; dimensione equ 16	rdx
; punto2 equ 12		rsi
; punto1 equ 8		rdi

creaMatriceDistanze:
		push rbp
		mov	rbp, rsp
		;pushaq     ;inizio

        mov eax, [rdi+d]
        cdq
        mov ecx, [rdi+m]
        idiv ecx
        mov rdx, rax

        mov r15, rdx
        sal r15, 2

        mov r14, [rdi+distanze_simmetriche]
        
        mov eax, [rdi+k]

        mov r8, rdi     ;input
        mov r9, rsi     ;codebook

        xor r13, r13
        xor r10, r10
        inc r10
forI:   cmp r10, rax
        jge endI
        mov rsi, r9

        xor r11, r11
forJ:   cmp r11, r10
        jge endJ

        mov edi, [r8+d]
        imul rdi, r10
        sal rdi, 2
        add rdi, r9

        xor r12, r12
        sub rcx, 4

forKQ:  cmp r12, rcx
        jge endKQ

        call distanza
        add rsi, r15
        add rdi, r15
        vmovss xmm11, xmm0

        call distanza
        add rsi, r15
        add rdi, r15
        vinsertps xmm11, xmm0, 00010000b
        
        call distanza
        add rsi, r15
        add rdi, r15
        vinsertps xmm11, xmm0, 00100000b
        
        call distanza
        add rsi, r15
        add rdi, r15
        vinsertps xmm11, xmm0, 00110000b
        
        vmovups [r14+4*r13], xmm11
        
        add r13, 4;count
        add r12, 4
        jmp forKQ
endKQ:  add rcx, 4
forKR:  cmp r12, rcx
        jge endKR

        call distanza
        vmovss [r14+4*r13], xmm0
        add rsi, r15
        add rdi, r15
        inc r13;count

        inc r12
        jmp forKR
endKR:   
        inc r11
        jmp forJ
endJ:   
        inc r10
        jmp forI
endI:   
        ;popaq		;fine
		mov	rsp, rbp
		pop	rbp
		ret