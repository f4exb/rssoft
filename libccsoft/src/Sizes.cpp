/*
     Copyright 2013 Edouard Griffiths <f4exb at free dot fr>

     This file is part of CCSoft. A Convolutional Codes Soft Decoding library

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin Street, Boston, MA  02110-1301  USA

     Objects sizes tester

*/

#include "CC_TreeNodeEdge.h"
#include "CC_TreeNodeEdge_FA.h"
#include "CC_TreeEdge.h"
#include "CC_TreeNode.h"

#include <iostream>


int main(int argc, char *argv[])
{
    std::cout << "unsigned int .....: " << sizeof(unsigned int) << " bytes" << std::endl;
    std::cout << "unsigned int* ....: " << sizeof(unsigned int*) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, empty>" << std::endl;
    std::cout << "TreeNode .........: " << sizeof(ccsoft::CC_TreeNode<unsigned int, unsigned int, ccsoft::CC_TreeEdgeTag_Empty>) << " bytes" << std::endl;
    std::cout << "TreeEdge .........: " << sizeof(ccsoft::CC_TreeEdge<unsigned int, unsigned int, ccsoft::CC_TreeEdgeTag_Empty>) << " bytes" << std::endl;
    std::cout << "TreeNodeEdge .....: " << sizeof(ccsoft::CC_TreeNodeEdge<unsigned int, unsigned int, ccsoft::CC_TreeNodeEdgeTag_Empty>) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, bool>" << std::endl;
    std::cout << "TreeNode .........: " << sizeof(ccsoft::CC_TreeNode<unsigned int, unsigned int, bool>) << " bytes" << std::endl;
    std::cout << "TreeEdge .........: " << sizeof(ccsoft::CC_TreeEdge<unsigned int, unsigned int, bool>) << " bytes" << std::endl;
    std::cout << "TreeNodeEdge .....: " << sizeof(ccsoft::CC_TreeNodeEdge<unsigned int, unsigned int, bool>) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, empty, 1>" << std::endl;
    std::cout << "TreeNodeEdge_FA ..: " << sizeof(ccsoft::CC_TreeNodeEdge_FA<unsigned int, unsigned int, ccsoft::CC_TreeNodeEdgeTag_Empty, 1>) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, bool, 1>" << std::endl;
    std::cout << "TreeNodeEdge_FA ..: " << sizeof(ccsoft::CC_TreeNodeEdge_FA<unsigned int, unsigned int, bool, 1>) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, bool, 2>" << std::endl;
    std::cout << "TreeNodeEdge_FA ..: " << sizeof(ccsoft::CC_TreeNodeEdge_FA<unsigned int, unsigned int, bool, 2>) << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "<unsigned int, unsigned int, bool, 3>" << std::endl;
    std::cout << "TreeNodeEdge_FA ..: " << sizeof(ccsoft::CC_TreeNodeEdge_FA<unsigned int, unsigned int, bool, 3>) << " bytes" << std::endl;

    return 0;
}
